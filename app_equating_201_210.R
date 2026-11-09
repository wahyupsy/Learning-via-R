# ================================================================
# TEST EQUATING 201-210 | R SHINY
# 4 pendekatan: Concurrent, Fixed Anchor, IRT Linking, CTT NEAT
# Install: install.packages(c('shiny','mirt','ggplot2','DT'))
# Jalankan: shiny::runApp('app_equating_201_210.R')
# ================================================================

library(shiny)
library(mirt)
library(ggplot2)
library(DT)

anchor_ids <- c('02','05','09','13','16','20','24','29','33','38','43',
                '47','48','50','52','66','67','68','73','83','84','85','86')
forms <- sprintf('form%d', 201:210)

# Kolom '02', '2', 'item02', 'I02', 'Q02', 'butir02' dikenali sebagai item 02.
# Kolom lain dianggap metadata, kecuali opsi "semua kolom respons" dipilih.
normalize_item <- function(x) {
  x <- trimws(tolower(gsub('^\\ufeff', '', x)))
  x <- gsub('[ ._-]', '', x)
  x <- sub('^(item|butir|soal|q|i)(?=[0-9]+$)', '', x, perl=TRUE)
  ifelse(grepl('^[0-9]+$', x), sprintf('I%03d', as.integer(x)), NA_character_)
}
anchor_names <- normalize_item(anchor_ids)

read_form <- function(path, form, sep=',', header=FALSE) {
  z <- read.csv(path, sep=sep, header=header, check.names=FALSE,
                stringsAsFactors=FALSE, na.strings=c('', 'NA', 'NaN', '.'))
  raw_names <- names(z)
  keys <- normalize_item(raw_names)
  chosen <- which(!is.na(keys))
  # CSV tanpa header: read.csv(header=FALSE) menghasilkan V1, V2, ...
  # Dalam mode ini nomor item diasumsikan sama dengan posisi kolom.
  # Dukung pula file yang telah dibaca ulang dengan header V1, V2, ...
  vcols <- grepl('^V[0-9]+$', raw_names)
  if (!length(chosen) && all(vcols)) {
    positions <- as.integer(sub('^V', '', raw_names))
    if (!identical(positions, seq_along(raw_names))) {
      stop(form, ': kolom V tidak berurutan dari V1. Tidak aman menentukan ID anchor.')
    }
    keys <- sprintf('I%03d', positions)
    chosen <- seq_along(raw_names)
  }
  if (!length(chosen)) stop(form, ': nama kolom tidak dikenali. Gunakan 01, item01, atau V1, V2, ...')
  if (anyDuplicated(keys[chosen])) stop(form, ': ID butir ganda setelah normalisasi. Periksa nama kolom.')
  missing_anchor <- setdiff(anchor_names, keys[chosen])
  if (length(missing_anchor)) stop(form, ': anchor tidak ditemukan: ',
                                   paste(sub('^I0*','',missing_anchor),collapse=', '),
                                   '. Untuk format V1..Vn, diperlukan sedikitnya 86 kolom, dan posisi kolom harus konsisten.')
  d <- z[, chosen, drop=FALSE]
  names(d) <- keys[chosen]
  for (nm in names(d)) {
    val <- trimws(as.character(d[[nm]]))
    val[val %in% c('', 'NA', 'NaN', '.')] <- NA
    bad <- !is.na(val) & !val %in% c('0','1')
    if (any(bad)) stop(form, ': kolom ', raw_names[chosen[match(nm, names(d))]],
                       ' memiliki nilai bukan 0/1. Lakukan rekoding terlebih dahulu.')
    d[[nm]] <- as.numeric(val)
  }
  # Nonanchor dengan ID sama di form berbeda merupakan item berbeda.
  names(d)[!names(d) %in% anchor_names] <- paste0(form, '_', names(d)[!names(d) %in% anchor_names])
  if (nrow(d) < 20) stop(form, ': minimal 20 peserta diperlukan.')
  if (any(rowSums(!is.na(d)) == 0)) stop(form, ': ada peserta tanpa respons sama sekali.')
  d
}

read_all <- function(folder, sep=',', header=FALSE) {
  paths <- file.path(folder, paste0(forms, '.csv'))
  if (any(!file.exists(paths))) stop('File tidak ditemukan: ', paste(basename(paths[!file.exists(paths)]),collapse=', '))
  d <- setNames(lapply(seq_along(forms), function(i) read_form(paths[i],forms[i],sep,header)), forms)
  for (f in forms) {
    if (any(vapply(d[[f]], function(x) all(is.na(x)) || length(unique(na.omit(x))) < 2, logical(1))))
      stop(f, ': ada butir tanpa variasi respons atau seluruhnya NA. Periksa data.')
  }
  d
}

item_type <- function(model) if (model=='Rasch') 'Rasch' else '2PL'
fit_single <- function(dat, model, pars=NULL) {
  mirt::mirt(dat, 1, itemtype=item_type(model), pars=pars,
             verbose=FALSE, technical=list(NCYCLES=500), SE=FALSE)
}
get_par <- function(fit) {
  co <- mirt::coef(fit, IRTpars=TRUE, simplify=TRUE)$items
  data.frame(item=rownames(co), a=as.numeric(co[,'a']),
             b=as.numeric(co[,'b']), row.names=NULL)
}

# Untuk perbandingan parameter, skala Rasch menggunakan A=1.
# MM: rata-rata diskriminasi dan kesukaran; MS: SD kesukaran.
# HB dan SL: minimisasi selisih kurva dengan bobot N(0,1).
prob <- function(theta, a, b) plogis(outer(theta,b,'-') *
                                    matrix(a,length(theta),length(a),byrow=TRUE))
link_irt <- function(px, py, method, model) {
  q <- merge(px[px$item %in% anchor_names,],
             py[py$item %in% anchor_names,], by='item', suffixes=c('.x','.y'))
  if(nrow(q)<3) stop('Minimal 3 anchor diperlukan untuk linking.')
  bx <- q$b.x; by <- q$b.y; ax <- q$a.x; ay <- q$a.y
  if(model=='Rasch') {
    A0 <- 1
  } else if(method=='MM') {
    A0 <- mean(ay)/mean(ax)
  } else {
    A0 <- sd(bx)/sd(by)
  }
  if (!is.finite(A0) || A0<=0) A0 <- 1
  B0 <- mean(bx)-A0*mean(by)
  if(method %in% c('MM','MS')) {
    if(method=='MS' && model!='Rasch') A0 <- sd(bx)/sd(by)
    if(method=='MS' && model=='Rasch') A0 <- 1
    return(c(A=A0, B=mean(bx)-A0*mean(by)))
  }
  th <- seq(-4,4,length.out=121)
  w <- dnorm(th); w <- w/sum(w)
  # Bandingkan pada skala X dan pada skala Y (kriteria simetris).
  loss <- function(v) {
    A <- if(model=='Rasch') 1 else exp(v[1])
    B <- if(model=='Rasch') v[1] else v[2]
    pxx <- prob(th,ax,bx)
    pyy_on_x <- prob((th-B)/A,ay,by)
    pyy <- prob(th,ay,by)
    pxx_on_y <- prob(A*th+B,ax,bx)
    if(method=='SL') {
      sum(w*(rowSums(pxx)-rowSums(pyy_on_x))^2) +
        sum(w*(rowSums(pyy)-rowSums(pxx_on_y))^2)
    } else {
      sum(w*rowSums((pxx-pyy_on_x)^2)) +
        sum(w*rowSums((pyy-pxx_on_y)^2))
    }
  }
  initial <- if(model=='Rasch') B0 else c(log(A0),B0)
  opt <- optim(initial,loss,method='BFGS',control=list(maxit=1500))
  A <- if(model=='Rasch') 1 else exp(opt$par[1])
  B <- if(model=='Rasch') opt$par[1] else opt$par[2]
  c(A=A,B=B)
}

transform_par <- function(par, A, B) {
  par$a <- par$a/A
  par$b <- A*par$b+B
  par
}

# CTT common-item chained linear equating, target Y -> baseline X.
# Y total -> Y anchor -> X anchor (identitas skor anchor) -> X total.
# Cocok sebagai demonstrasi NEAT; perlu asumsi invariansi anchor.
ctt_equate <- function(dx, dy) {
  sx <- rowSums(dx, na.rm=FALSE); sy <- rowSums(dy, na.rm=FALSE)
  ax <- rowSums(dx[,anchor_names,drop=FALSE],na.rm=FALSE)
  ay <- rowSums(dy[,anchor_names,drop=FALSE],na.rm=FALSE)
  # Complete cases untuk CTT agar denominator skor total konsisten.
  okx <- is.finite(sx) & is.finite(ax)
  oky <- is.finite(sy) & is.finite(ay)
  if(sum(okx)<20 || sum(oky)<20) stop('CTT: minimal 20 peserta lengkap per form.')
  sx <- sx[okx]; sy <- sy[oky]; ax <- ax[okx]; ay <- ay[oky]
  if(sd(ax)==0 || sd(ay)==0 || sd(sy)==0) stop('CTT: SD skor nol.')
  y_grid <- 0:ncol(dy)
  y_to_ay <- mean(ay) + sd(ay)/sd(sy)*(y_grid-mean(sy))
  equated <- mean(sx) + sd(sx)/sd(ax)*(y_to_ay-mean(ax))
  data.frame(skor_Y=y_grid, skor_setara_X=equated,
             skor_setara_X_dibatasi=pmax(0,pmin(ncol(dx),equated)))
}

# Concurrent multiple-group calibration: anchor parameters invariant,
# mean/variance kelompok target bebas, baseline X mean=0, SD=1.
fit_concurrent <- function(d, model) {
  all_items <- unique(unlist(lapply(d,names)))
  blocks <- lapply(d,function(x) {
    z <- matrix(NA_real_,nrow(x),length(all_items),dimnames=list(NULL,all_items))
    z[,names(x)] <- as.matrix(x)
    z
  })
  stacked <- do.call(rbind,blocks)
  grp <- factor(rep(names(d),lengths(lapply(d,function(x) seq_len(nrow(x))))), levels=forms)
  fit <- mirt::multipleGroup(stacked, 1, group=grp,
            itemtype=item_type(model),
            invariance=c(anchor_names,'free_means','free_var'),
            verbose=FALSE, technical=list(NCYCLES=700), SE=FALSE)
  out <- lapply(forms,function(f) {
    cf <- mirt::coef(fit,group=f,IRTpars=TRUE,simplify=TRUE)$items
    p <- data.frame(item=rownames(cf),a=as.numeric(cf[,'a']),b=as.numeric(cf[,'b']))
    p[p$item %in% names(d[[f]]),,drop=FALSE]
  })
  names(out) <- forms
  list(fit=fit, pars=out)
}

# Fixed anchor: mirt baseline -> mod2values target; d dan a1 anchor fixed.
fit_fixed <- function(d, model, baseline_fit) {
  vbase <- mirt::mod2values(baseline_fit)
  out <- list(form201=get_par(baseline_fit))
  for(f in setdiff(forms,'form201')) {
    vt <- mirt::mirt(d[[f]],1,itemtype=item_type(model),pars='values',verbose=FALSE)
    for(id in anchor_names) for(pn in c('a1','d')) {
      src <- vbase[vbase$item==id & vbase$name==pn, 'value']
      idx <- which(vt$item==id & vt$name==pn)
      if(length(src)==1 && length(idx)==1) {
        vt$value[idx] <- src
        vt$est[idx] <- FALSE
      }
    }
    fit <- fit_single(d[[f]],model,pars=vt)
    out[[f]] <- get_par(fit)
  }
  out
}

ui <- fluidPage(
  titlePanel('Test Equating 201–210 | Baseline form201'),
  sidebarLayout(
    sidebarPanel(width=3,
      textInput('folder','Folder dataset (Windows)',
        value='D:/OneDrive/Github/Learning via R/dataset equating'),
      selectInput('sep','Pemisah CSV',choices=c('Koma'=',' ,'Titik koma'=';'),selected=','),
      checkboxInput('has_header','CSV memiliki baris header',value=FALSE),
      selectInput('model','Model IRT',choices=c('Rasch','2PL'),selected='Rasch'),
      selectInput('approach','Pendekatan',choices=c(
        'Concurrent calibration'='concurrent',
        'Fixed item calibration'='fixed',
        'Transformasi linier IRT'='link',
        'Classical Test Equating (CTT)'='ctt')),
      conditionalPanel("input.approach == 'link'",
        selectInput('method','Metode linking untuk visualisasi',
                    choices=c('MM','MS','SL','HB'),selected='SL')),
      selectInput('target','Form target',choices=setNames(forms[-1],forms[-1])),
      actionButton('run','Jalankan equating',class='btn-primary'),
      hr(), strong('Anchor (23 butir)'),
      p(paste(anchor_ids,collapse=', ')),
      p('Untuk CSV tanpa header: V1 = item 01, V2 = item 02, ..., V86 = item 86. Pastikan posisi anchor sama di semua form.'),
      p('Butir nonanchor bernomor sama pada form berbeda diperlakukan sebagai butir berbeda.'),
      downloadButton('download','Unduh hasil CSV')),
    mainPanel(width=9,
      tabsetPanel(
        tabPanel('1. Data & skor',plotOutput('rawplot',height=380),DTOutput('data_summary')),
        tabPanel('2. Parameter butir',DTOutput('params')),
        tabPanel('3. Perubahan parameter',plotOutput('parplot',height=430)),
        tabPanel('4. ICC + IIF anchor',selectInput('anchor','Butir anchor',choices=anchor_names),plotOutput('anchorplot',height=470)),
        tabPanel('5. ICC + IIF semua',uiOutput('item_ui'),plotOutput('itemplot',height=470)),
        tabPanel('6. Skor CTT',plotOutput('cttplot',height=410),DTOutput('ctttable')),
        tabPanel('7. TCC',plotOutput('tccplot',height=440)),
        tabPanel('8. MM MS SL HB',DTOutput('links'),plotOutput('linkplot',height=420)),
        tabPanel('9. Semua form',DTOutput('allforms')),
        tabPanel('10. Petunjuk',
          h4('Desain penyetaraan'),
          p('Form201 adalah baseline. Anchor identik antarform. Nonanchor diberi nama unik per form.'),
          p('Concurrent: seluruh form dikalibrasi bersama; anchor invariant; mean dan varians kelompok target bebas.'),
          p('Fixed: parameter anchor pada baseline dikunci dalam kalibrasi masing-masing target.'),
          p('Linking: kalibrasi terpisah, kemudian MM/MS/SL/HB pada anchor. A dan B mengubah theta target ke theta baseline.'),
          p('CTT: chained linear equating dengan anchor, menggunakan peserta dengan respons lengkap.'),
          p('Grafik ICC/IIF dan TCC membandingkan parameter baseline dan target pada skala baseline.'),
          p('Pada Rasch A=1; metode MS dan MM sama-sama menjadi mean shift, sementara SL/HB dapat berbeda.'),
          p('Tidak semua form harus memiliki jumlah butir yang sama. Interpretasi skor total yang setara perlu mempertimbangkan panjang tes.'),
          p('Keterbatasan: hasil tidak otomatis membuktikan invariansi anchor, validitas desain NEAT, atau ketepatan model.'))
      )
    )
  )
)

server <- function(input,output,session) {
  dataset <- eventReactive(input$run, {
    withProgress(message='Membaca dataset...',value=0.1, {
      d <- read_all(input$folder,input$sep,input$has_header)
      incProgress(0.2)
      d
    })
  })
  result <- eventReactive(input$run, {
    d <- dataset(); model <- isolate(input$model); approach <- isolate(input$approach)
    withProgress(message='Kalibrasi / equating sedang diproses...',value=0.1, {
      if(approach=='ctt') {
        maps <- lapply(forms[-1],function(f) ctt_equate(d$form201,d[[f]]))
        names(maps) <- forms[-1]
        incProgress(0.8)
        return(list(d=d,approach=approach,maps=maps,par=NULL,link=NULL))
      }
      if(approach=='concurrent') {
        par <- fit_concurrent(d,model)$pars
        incProgress(0.8)
        return(list(d=d,approach=approach,par=par,link=NULL,maps=NULL))
      }
      base <- fit_single(d$form201,model)
      incProgress(0.15)
      if(approach=='fixed') {
        par <- fit_fixed(d,model,base)
        incProgress(0.7)
        return(list(d=d,approach=approach,par=par,link=NULL,maps=NULL))
      }
      raw <- lapply(d,function(x) get_par(fit_single(x,model)))
      incProgress(0.4)
      links <- do.call(rbind,lapply(forms[-1],function(f) {
        do.call(rbind,lapply(c('MM','MS','SL','HB'),function(m) {
          ab <- link_irt(raw$form201,raw[[f]],m,model)
          data.frame(form=f,method=m,A=unname(ab['A']),B=unname(ab['B']))
        }))
      }))
      incProgress(0.4)
      list(d=d,approach=approach,raw=raw,link=links,par=NULL,maps=NULL)
    })
  })
  output$item_ui <- renderUI({
    r <- result(); req(r)
    selectInput('item','Butir pada form target',choices=names(r$d[[input$target]]))
  })
  get_pair <- reactive({
    r <- result(); req(r)
    if(r$approach=='ctt') return(NULL)
    x <- if(r$approach=='link') r$raw$form201 else r$par$form201
    y <- if(r$approach=='link') {
      ab <- subset(r$link,form==input$target & method==input$method)
      transform_par(r$raw[[input$target]],ab$A[1],ab$B[1])
    } else r$par[[input$target]]
    list(x=x,y=y)
  })
  output$rawplot <- renderPlot({
    d <- result()$d
    v <- do.call(rbind,lapply(names(d),function(f)
      data.frame(form=f,skor=rowSums(d[[f]],na.rm=FALSE))))
    ggplot(v,aes(skor,fill=form))+geom_histogram(bins=30,alpha=.7)+
      facet_wrap(~form,scales='free')+theme_minimal()+
      labs(x='Skor mentah',y='Frekuensi',title='Distribusi skor mentah per form')+
      theme(legend.position='none')
  })
  output$data_summary <- renderDT({
    d <- result()$d
    datatable(do.call(rbind,lapply(names(d),function(f) data.frame(
      form=f,n=nrow(d[[f]]),n_item=ncol(d[[f]]),n_anchor=length(anchor_names),
      rerata_skor=mean(rowSums(d[[f]],na.rm=FALSE),na.rm=TRUE),
      sd_skor=sd(rowSums(d[[f]],na.rm=FALSE),na.rm=TRUE)))),options=list(pageLength=10))
  })
  output$params <- renderDT({
    p <- get_pair(); validate(need(!is.null(p),'CTT tidak mengestimasi parameter IRT.'))
    x <- transform(p$x,form='form201'); y <- transform(p$y,form=input$target)
    z <- rbind(x,y); z$anchor <- z$item %in% anchor_names
    datatable(z,options=list(pageLength=25,scrollX=TRUE))
  })
  output$parplot <- renderPlot({
    p <- get_pair(); validate(need(!is.null(p),'Pilih pendekatan IRT untuk grafik parameter.'))
    q <- merge(p$x[p$x$item %in% anchor_names,],p$y[p$y$item %in% anchor_names,],
               by='item',suffixes=c('_X','_Y'))
    ggplot(q,aes(b_X,b_Y,label=item))+geom_point(size=2,color='#2673b8')+
      geom_abline(slope=1,intercept=0,linetype=2)+
      geom_text(vjust=-.7,size=3,check_overlap=TRUE)+coord_equal()+theme_minimal()+
      labs(x='b baseline form201',y=paste('b',input$target,'pada skala baseline'),
           title='Kesulitan anchor: garis putus-putus = kesetaraan sempurna')
  })
  curve_plot <- function(item) {
    p <- get_pair(); validate(need(!is.null(p),'ICC/IIF tersedia untuk IRT.'))
    th <- seq(-4,4,length.out=161)
    rows <- list()
    for(f in c('form201',input$target)) {
      par <- if(f=='form201') p$x else p$y
      z <- par[par$item==item,,drop=FALSE]
      if(nrow(z)==0) next
      pp <- plogis(z$a[1]*(th-z$b[1]))
      rows[[f]] <- data.frame(theta=th,ICC=pp,IIF=z$a[1]^2*pp*(1-pp),form=f)
    }
    validate(need(length(rows)>0,'Butir tidak ditemukan pada form yang dipilih.'))
    v <- do.call(rbind,rows)
    # IIF diskalakan hanya untuk overlay visual; sumbu y kedua menunjukkan nilai asli.
    k <- max(v$IIF,na.rm=TRUE)
    v$IIFscaled <- v$IIF/k
    ggplot(v,aes(theta,color=form))+geom_line(aes(y=ICC,linetype='ICC'),linewidth=1)+
      geom_line(aes(y=IIFscaled,linetype='IIF (dinormalisasi)'),linewidth=.9)+
      scale_y_continuous('Peluang benar / IIF relatif',sec.axis=sec_axis(~.*k,name='IIF asli'))+
      scale_linetype_manual(values=c('ICC'='solid','IIF (dinormalisasi)'='dashed'))+
      theme_minimal()+labs(title=paste('ICC dan IIF:',item),x='Theta pada skala form201',linetype='Kurva')
  }
  output$anchorplot <- renderPlot(curve_plot(input$anchor))
  output$itemplot <- renderPlot(curve_plot(input$item))
  output$cttplot <- renderPlot({
    r <- result(); validate(need(r$approach=='ctt','Tab ini khusus CTT.'))
    m <- r$maps[[input$target]]
    ggplot(m,aes(skor_Y,skor_setara_X))+geom_line(linewidth=1,color='#2673b8')+
      geom_point(size=1)+theme_minimal()+labs(x='Skor mentah target',
      y='Skor setara pada baseline',title=paste('CTT chained linear:',input$target,'ke form201'))
  })
  output$ctttable <- renderDT({
    r <- result(); validate(need(r$approach=='ctt','Pilih pendekatan CTT.'))
    datatable(r$maps[[input$target]],options=list(pageLength=15))
  })
  output$tccplot <- renderPlot({
    p <- get_pair(); validate(need(!is.null(p),'TCC tersedia untuk IRT.'))
    th <- seq(-4,4,length.out=161)
    v <- do.call(rbind,lapply(c('form201',input$target),function(f) {
      z <- if(f=='form201') p$x else p$y
      data.frame(theta=th,TCC=rowSums(prob(th,z$a,z$b)),form=f)
    }))
    ggplot(v,aes(theta,TCC,color=form))+geom_line(linewidth=1)+theme_minimal()+
      labs(x='Theta skala baseline',y='Ekspektasi skor mentah',title='Test Characteristic Curve')
  })
  output$links <- renderDT({
    r <- result(); validate(need(r$approach=='link','Pilih Transformasi Linier IRT.'))
    datatable(subset(r$link,form==input$target),options=list(dom='t'),rownames=FALSE)
  })
  output$linkplot <- renderPlot({
    r <- result(); validate(need(r$approach=='link','Pilih Transformasi Linier IRT.'))
    q <- subset(r$link,form==input$target)
    ggplot(q,aes(method,B,fill=method))+geom_col(width=.6)+theme_minimal()+
      labs(x='Metode linking',y='Konstanta B',title=paste('Perbandingan konstanta B:',input$target))+
      theme(legend.position='none')
  })
  output$allforms <- renderDT({
    r <- result()
    if(r$approach=='link') return(datatable(r$link,options=list(pageLength=36)))
    if(r$approach=='ctt') {
      z <- do.call(rbind,lapply(names(r$maps),function(f) {
        m <- r$maps[[f]]
        data.frame(form=f,skor_min_setara=m$skor_setara_X[1],
                   skor_max_setara=tail(m$skor_setara_X,1))
      }))
      return(datatable(z,options=list(pageLength=10)))
    }
    z <- do.call(rbind,lapply(names(r$par),function(f) {
      p <- r$par[[f]]
      data.frame(form=f,jumlah_item=nrow(p),rerata_a=mean(p$a),rerata_b=mean(p$b))
    }))
    datatable(z,options=list(pageLength=10))
  })
  output$download <- downloadHandler(
    filename=function() paste0('equating_201_210_',Sys.Date(),'.csv'),
    content=function(file) {
      r <- result()
      if(r$approach=='link') z <- r$link
      else if(r$approach=='ctt') z <- do.call(rbind,lapply(names(r$maps),function(f)
        cbind(form=f,r$maps[[f]])))
      else z <- do.call(rbind,lapply(names(r$par),function(f) cbind(form=f,r$par[[f]])))
      write.csv(z,file,row.names=FALSE)
    }
  )
  observeEvent(result(),{
    showNotification('Analisis selesai. Semua form telah diproses.',type='message',duration=5)
  })
}

shinyApp(ui,server)
