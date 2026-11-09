# ============================================================
# TEST EQUATING INTERAKTIF - 10 FORM, 4 PENDEKATAN
# Sumber: adaptasi skrip equating berpasangan yang telah berhasil
# Jalankan: shiny::runApp('app_equating_interaktif_TAM_Rasch_filter_baris_NA.R')
# ============================================================
required <- c('shiny','TAM','ggplot2','DT')
missing <- required[!vapply(required, requireNamespace, logical(1), quietly=TRUE)]
if(length(missing)) stop('Instal paket: install.packages(c(',paste(sprintf('"%s"',missing),collapse=','),'))')
library(shiny); library(TAM); library(ggplot2); library(DT)

folder_default <- 'D:/OneDrive/Github/Learning via R/dataset equating'
form_ids <- sprintf('form%d',201:210)
anchor <- sprintf('%02d',c(2,5,9,13,16,20,24,29,33,38,43,47,48,50,52,66,67,68,73,83,84,85,86))
subtests <- list('Subtes 1 (1–30)' = 1:30,
                 'Subtes 2 (31–60)' = 31:60,
                 'Subtes 3 (61–90)' = 61:90)
subtest_options <- setNames(names(subtests), names(subtests))
anchor_for_subtest <- function(items) intersect(anchor, sprintf('%02d',items))

read_form <- function(path, subtest_items) {
  if(!file.exists(path)) stop('File tidak ditemukan: ',path)
  z <- read.csv(path,header=FALSE,check.names=FALSE,na.strings=c('','NA'))
  if (ncol(z) < max(subtest_items)) {
    stop(basename(path), ' hanya memiliki ', ncol(z),
         ' kolom; subtes ini memerlukan kolom sampai ', max(subtest_items), '.')
  }
  names(z) <- sprintf('%02d',seq_len(ncol(z)))
  z <- z[, sprintf('%02d',subtest_items), drop=FALSE]
  z[] <- lapply(z,function(v) suppressWarnings(as.numeric(as.character(v))))
  v <- unlist(z,use.names=FALSE)
  if(any(!is.na(v)& !v %in% c(0,1))) stop('Respons bukan 0/1: ',basename(path))
  # Hapus hanya peserta dengan semua respons NA; jangan ubah NA sebagian menjadi nol.
  empty_rows <- rowSums(!is.na(z)) == 0L
  removed <- sum(empty_rows)
  if (removed > 0L) {
    message('[', basename(path), '] Menghapus ', removed, ' baris seluruhnya NA dari ', nrow(z), ' baris.')
    z <- z[!empty_rows, , drop=FALSE]
  }
  if (nrow(z) < 2L) stop('Terlalu sedikit peserta valid setelah pembersihan: ', basename(path))
  attr(z, 'removed_all_na') <- removed
  attr(z, 'original_n') <- nrow(z) + removed
  if(any(vapply(z,function(v) length(unique(v[!is.na(v)]))<2,logical(1))))
    stop('Ada butir dengan respons konstan pada ',basename(path),' (estimasi memerlukan variasi).')
  z
}
# Sampel 5% sampai 80% dari peserta valid tiap form.
# Seed tetap agar sampel konsisten saat aplikasi dijalankan ulang.
sample_form <- function(d, percentage, seed) {
  n_total <- nrow(d)
  n_keep <- max(2L, min(n_total, as.integer(round(n_total * percentage / 100))))
  set.seed(seed)
  selected <- sort(sample.int(n_total, size=n_keep, replace=FALSE))
  out <- d[selected, , drop=FALSE]
  attr(out, 'removed_all_na') <- attr(d, 'removed_all_na')
  attr(out, 'original_n') <- attr(d, 'original_n')
  attr(out, 'valid_n') <- n_total
  attr(out, 'sample_pct') <- percentage
  # Hanya validasi: jika subsampling menghasilkan butir tanpa variasi,
  # berikan pesan yang dapat ditindaklanjuti.
  constant <- names(out)[vapply(out, function(z) length(unique(z[!is.na(z)])) < 2L, logical(1))]
  if (length(constant)) {
    stop('Subsampling ',percentage,'% menghasilkan butir konstan pada ',
         paste(head(constant,12),collapse=', '),
         '. Coba naikkan persentase sampel (misalnya 10% atau 20%).')
  }
  out
}

# Wrapper: nama butir aman untuk TAM, tetapi nomor butir asli tetap disimpan.
tam_fit <- function(d, ...) {
  ids <- colnames(d)
  names(d) <- sprintf('Item%03d', seq_along(ids))
  fit <- TAM::tam.mml(resp=d, verbose=FALSE, ...)
  fit$equating_item_ids <- ids
  fit
}

# Ekstraksi kesukaran Rasch dari TAM. Parameter xsi adalah kesukaran b.
get_ab <- function(fit) {
  ids <- fit$equating_item_ids
  if (is.null(ids)) stop('Metadata nama butir hasil TAM tidak ditemukan.')
  # Untuk Rasch dikotomis, urutan parameter xsi mengikuti urutan butir.
  # Parameter xsi pada TAM merupakan parameter kesukaran b.
  b <- as.numeric(fit$xsi$xsi)
  if (length(b) != length(ids)) {
    stop('Jumlah parameter TAM (', length(b), ') tidak sama dengan jumlah butir (',
         length(ids), '). Periksa spesifikasi model Rasch dikotomis.')
  }
  if (any(!is.finite(b))) stop('Parameter b TAM mengandung nilai tidak terhingga / NA.')
  data.frame(item=ids, a=rep(1, length(b)), b=b, row.names=NULL)
}

# Bulatkan kolom numerik di SEMUA tabel hasil (3 digit).
round3 <- function(d) {
  d[] <- lapply(d, function(x) if (is.numeric(x)) round(x, 3) else x)
  d
}

show_table <- function(d, pageLength=12, dom=NULL) {
  ops <- list(pageLength=pageLength, scrollX=TRUE)
  if (!is.null(dom)) ops$dom <- dom
  tab <- DT::datatable(round3(d), options=ops, rownames=FALSE)
  cols <- names(d)[vapply(d, is.numeric, logical(1))]
  if(length(cols)) tab <- DT::formatRound(tab, columns=cols, digits=3)
  tab
}

prob <- function(theta,par) {
  plogis(outer(theta,par$b,'-') * matrix(par$a,length(theta),nrow(par),byrow=TRUE))
}
info <- function(theta,par) {
  p <- prob(theta,par)
  p*(1-p)*matrix(par$a^2,length(theta),nrow(par),byrow=TRUE)
}
# Linking lokasi Rasch: Y -> X, A selalu 1.
# MM dan MS akan sama ketika slope sengaja dibatasi 1.
link_ab <- function(px,py,method,anchors,grid=seq(-4,4,length.out=121)) {
  ix <- match(anchors,px$item); iy <- match(anchors,py$item)
  if (anyNA(ix) || anyNA(iy)) stop('Anchor tidak lengkap pada parameter TAM.')
  bx <- px$b[ix]; by <- py$b[iy]
  if (method %in% c('MM','MS')) {
    B <- mean(bx)-mean(by)
  } else {
    target <- plogis(outer(grid,bx,'-'))
    objective <- function(B) {
      candidate <- plogis(outer(grid,by+B,'-'))
      if (method == 'SL') {
        sum((rowSums(target)-rowSums(candidate))^2)
      } else {
        sum((target-candidate)^2)
      }
    }
    opt <- optimize(objective, interval=c(-8,8))
    B <- opt$minimum
  }
  data.frame(method=method,A=1,B=B)
}

# Parameter anchor dari X dikunci sewaktu kalibrasi Y.
fixed_calibrate <- function(Xfit,Y,anchors) {
  pX <- get_ab(Xfit)
  idx <- match(anchors,names(Y))
  val <- pX$b[match(anchors,pX$item)]
  if (anyNA(idx) || anyNA(val)) stop('Anchor tidak lengkap untuk fixed calibration.')
  tam_fit(Y, xsi.fixed=cbind(idx,val))
}

ctt_equate <- function(X,Y,anchors) {
  # CTT menggunakan kasus lengkap per form (tidak menganggap NA sebagai salah).
  X <- X[complete.cases(X),,drop=FALSE]
  Y <- Y[complete.cases(Y),,drop=FALSE]
  if(nrow(X)<10 || nrow(Y)<10) stop('CTT: kurang dari 10 respons lengkap per form; diperlukan metode missing-data khusus.')
  tx <- rowSums(X); ty <- rowSums(Y)
  ax <- rowSums(X[,anchors,drop=FALSE]); ay <- rowSums(Y[,anchors,drop=FALSE])
  if(any(!is.finite(c(sd(tx),sd(ty),sd(ax),sd(ay)))) || any(c(sd(tx),sd(ty),sd(ax),sd(ay))<=0)) stop('SD skor nol; CTT linear tidak dapat dihitung.')
  slope <- sd(tx)/sd(ax)*sd(ay)/sd(ty)
  intercept <- mean(tx)+sd(tx)/sd(ax)*(mean(ay)-mean(ax))-slope*mean(ty)
  data.frame(score_Y=0:ncol(Y),score_X_equated=intercept+slope*(0:ncol(Y)))
}

# Skor prorata hanya untuk visualisasi deskriptif; bukan skor CTT equating.
score_display <- function(z) {
  n <- rowSums(!is.na(z)); k <- ncol(z)
  out <- rowSums(z,na.rm=TRUE)*k/pmax(n,1)
  out[n==0] <- NA_real_
  out
}

run_pair <- function(X,Y,anchors) {
  if (length(anchors) < 2L) stop('Minimal dua anchor diperlukan untuk equating.')
  # 1. Kalibrasi Rasch terpisah dengan TAM
  fitX <- tam_fit(X)
  fitY <- tam_fit(Y)
  parX <- get_ab(fitX); parY <- get_ab(fitY)

  # 2. Concurrent: parameter anchor benar-benar dibagi, non-anchor
  # dipisahkan sebagai item yang berbeda; NA adalah planned missing.
  nonX <- setdiff(names(X),anchors); nonY <- setdiff(names(Y),anchors)
  Xc <- X[,c(anchors,nonX),drop=FALSE]; names(Xc) <- c(anchors,paste0('X_',nonX))
  Yc <- Y[,c(anchors,nonY),drop=FALSE]; names(Yc) <- c(anchors,paste0('Y_',nonY))
  common <- union(names(Xc),names(Yc))
  pad <- function(d) {
    z <- matrix(NA_real_,nrow(d),length(common),dimnames=list(NULL,common))
    z[,names(d)] <- as.matrix(d)
    as.data.frame(z,check.names=FALSE)
  }
  joined <- rbind(pad(Xc),pad(Yc))
  grp <- c(rep(1,nrow(X)),rep(2,nrow(Y)))
  empty_joined <- rowSums(!is.na(joined)) == 0L
  if (any(empty_joined)) {
    message('[Concurrent] Menghapus ', sum(empty_joined), ' baris seluruhnya NA pada matriks gabungan.')
    joined <- joined[!empty_joined, , drop=FALSE]
    grp <- grp[!empty_joined]
  }
  if (!all(table(factor(grp, levels=1:2)) >= 2L))
    stop('Jumlah peserta pada salah satu form terlalu sedikit setelah membersihkan baris seluruhnya NA.')
  concurrent <- tryCatch({
    fit <- tam_fit(joined, group=grp)
    pars <- get_ab(fit)
    cx <- pars[pars$item %in% names(Xc),,drop=FALSE]
    cy <- pars[pars$item %in% names(Yc),,drop=FALSE]
    cx$item <- sub('^X_','',cx$item)
    cy$item <- sub('^Y_','',cy$item)
    list(X=cx,Y=cy,error=NULL)
  },error=function(e) list(X=NULL,Y=NULL,error=conditionMessage(e)))

  # 3. Fixed anchor
  fixed <- tryCatch(list(par=get_ab(fixed_calibrate(fitX,Y,anchors)),error=NULL),
                    error=function(e) list(par=NULL,error=conditionMessage(e)))

  # 4. Transformasi linear, slope A=1 pada Rasch
  links <- do.call(rbind,lapply(c('MM','MS','SL','HB'),function(m) link_ab(parX,parY,m,anchors)))
  transformed <- lapply(seq_len(nrow(links)),function(j) {
    z <- parY
    z$b <- z$b+links$B[j]
    z
  })
  names(transformed) <- links$method

  # 5. CTT equating (skor total)
  ctt <- tryCatch(list(tab=ctt_equate(X,Y,anchors),error=NULL),
                  error=function(e) list(tab=NULL,error=conditionMessage(e)))
  list(X=X,Y=Y,anchors=anchors,parX=parX,parY=parY,concurrent=concurrent,fixed=fixed,
       links=links,transformed=transformed,ctt=ctt)
}

ui <- fluidPage(
  titlePanel('Laboratorium Test Equating Rasch (TAM) — Dua Form per Subtes'),
  sidebarLayout(
    sidebarPanel(width=3,
      textInput('folder','Folder CSV',value=folder_default),
      selectInput('formX','Form acuan (X)',choices=form_ids,selected='form201'),
      selectInput('formY','Form target (Y)',choices=form_ids,selected='form202'),
      selectInput('subtest','Pilih subtes',choices=subtest_options,
                  selected='Subtes 1 (1–30)'),
      helpText('Setiap subtes dikalibrasi dan disetarakan secara terpisah; anchor dibatasi pada subtes terpilih.'),
      sliderInput('sample_pct','Persentase sampel dipakai per form',
                  min=5,max=80,value=5,step=5,post='%' ),
      helpText('Subsampling acak tanpa pengembalian, seed tetap 2026. Klik Jalankan equating untuk menerapkan ukuran sampel.'),
      helpText('Model IRT: Rasch (TAM::tam.mml)'),
      selectInput('approach','Pendekatan yang divisualisasikan',
                  choices=c('Transformasi linier'='Linear',
                            'Fixed item calibration'='Fixed',
                            'Concurrent calibration'='Concurrent',
                            'CTT linear'='CTT'),selected='Linear'),
      conditionalPanel("input.approach == 'Linear'",
                       selectInput('method','Metode linking',choices=c('MM','MS','SL','HB'),selected='SL')),
      radioButtons('icc_items','Set item untuk plot ICC',choices=c('Anchor subtes'='anchor','Semua item dalam subtes'='common'),selected='anchor'),
      checkboxInput('anchor_only','Tampilkan hanya anchor pada grafik parameter',FALSE),
      actionButton('run','Jalankan equating',class='btn-primary'),
      hr(),helpText('Klik Jalankan equating setelah mengganti pasangan form. Pilihan pendekatan, metode, dan item dapat diganti tanpa kalibrasi ulang.'),
      verbatimTextOutput('status')
    ),
    mainPanel(width=9,
      tabsetPanel(
        tabPanel('1. Statistik item',
                 h4('Statistik item anchor — sebelum equating'),DTOutput('statAnchor'),
                 h4('Statistik seluruh item — sebelum equating'),DTOutput('statAll'),
                 h4('Statistik item anchor — setelah equating'),DTOutput('statAfter'),
                 p('Delta = rerata kesukaran Y dikurangi X; r_par = korelasi Pearson parameter b anchor. Pada hasil setelah equating, delta menggunakan Y setelah dan X baseline.')),
        tabPanel('2. Distribusi skor',plotOutput('rawPlot',height='620px'),DTOutput('rawTable')),
        tabPanel('3. Parameter butir',
                 p('Form X tetap sebagai baseline. Hanya form Y yang memiliki kolom sebelum dan sesudah equating. Untuk CTT tidak tersedia parameter item setelah equating.'),
                 DTOutput('paramTable')),
        tabPanel('4. Sebelum vs sesudah',plotOutput('paramPlot',height='450px'),DTOutput('diffTable')),
        tabPanel('5. ICC sebelum equating',plotOutput('iccPrePlot',height='1100px')),
        tabPanel('6. ICC setelah equating',plotOutput('iccPostPlot',height='1100px')),
        tabPanel('7. Distribusi skor setara',plotOutput('equatedPlot',height='450px'),DTOutput('conversionTable')),
        tabPanel('8. TCC',h4('TCC sebelum equating'),plotOutput('tccBefore',height='400px'),h4('TCC setelah equating'),plotOutput('tccAfter',height='400px')),
        tabPanel('9. MM MS SL HB',DTOutput('linkTable'),plotOutput('linkPlot',height='410px')),
        tabPanel('10. Ringkasan',DTOutput('summaryTable'),verbatimTextOutput('notes'),DTOutput('missingTable')), 
        tabPanel('11. Panduan',
          h4('Cara membaca grafik'),
          p('Kurva X adalah acuan. Kurva Y sebelum linking masih berada pada skala hasil kalibrasi Y; perbandingan langsung harus ditafsirkan hati-hati. Kurva Y setelah linking berada pada metrik X.'),
          p('Tab 5 menampilkan ICC sebelum equating. Tab 6 menampilkan ICC setelah equating. Tab 8 menyajikan TCC sebelum di atas dan sesudah di bawah.'),
          p('Semua tabel dan grafik mengikuti subtes terpilih (1–30, 31–60, atau 61–90). Anchor dibatasi pada subtes tersebut. Pilihan subtes baru berlaku setelah klik Jalankan equating.'),
          p('Set item untuk ICC dapat dipilih: hanya anchor atau semua nomor item pada subtes. Format grafik menggunakan facet per item.'),
          p('Concurrent: anchor diperlakukan invariant dalam kalibrasi bersama. Fixed: parameter anchor X dikunci saat kalibrasi Y. Linear: Y dikalibrasi sendiri lalu ditransformasi ke X.'),
          p('CTT: chained linear equating dengan anchor sebagai penghubung. Ini penyetaraan skor total, bukan transformasi parameter IRT.'),
          p('Perbandingan skor mentah sebelum/sesudah hanya bermakna setelah dibuat fungsi konversi. Untuk IRT digunakan expected-score true-score equating melalui TCC; CTT menggunakan konversi linear.'),
          p('Asumsi: posisi kolom anchor sama di seluruh form, CSV tanpa header, respons 0/1/NA. Hanya Rasch dengan TAM; A=1. MM dan MS sama karena linking lokasi saja. Skor deskriptif diprorata jika ada NA; CTT memakai peserta dengan respons lengkap saja.')
        )
      )
    )
  )
)

server <- function(input,output,session) {
  result <- eventReactive(input$run,{
    if(identical(input$formX,input$formY)) stop('Form X dan Y harus berbeda.')
    withProgress(message='Kalibrasi dan equating',value=0,{
      active_items <- subtests[[input$subtest]]
      if (is.null(active_items)) stop('Subtes tidak dikenali.')
      anchors <- anchor_for_subtest(active_items)
      X_full <- read_form(file.path(input$folder,paste0(input$formX,'.csv')), active_items)
      Y_full <- read_form(file.path(input$folder,paste0(input$formY,'.csv')), active_items)
      X <- sample_form(X_full, input$sample_pct, seed=2026L)
      Y <- sample_form(Y_full, input$sample_pct, seed=2027L)
      incProgress(.15,detail='Data dibaca dan sampel dipilih')
      ans <- run_pair(X,Y,anchors)
      ans$subtest <- input$subtest
      ans$formX <- input$formX
      ans$formY <- input$formY
      ans$sample_pct <- input$sample_pct
      incProgress(.85,detail='Selesai')
      ans
    })
  },ignoreInit=TRUE)
  output$status <- renderText({
    if(input$run==0) return('Pilih pasangan form, lalu klik Jalankan equating.')
    r <- tryCatch(result(),error=function(e)e)
    if(inherits(r,'error')) return(paste('ERROR:',conditionMessage(r)))
    paste0('X: ',r$formX,' (N valid=',attr(r$X,'valid_n'),', N sampel=',nrow(r$X),
           ', ',ncol(r$X),' butir; baris NA semua dibuang=',attr(r$X,'removed_all_na'),')\n',
           'Y: ',r$formY,' (N valid=',attr(r$Y,'valid_n'),', N sampel=',nrow(r$Y),
           ', ',ncol(r$Y),' butir; baris NA semua dibuang=',attr(r$Y,'removed_all_na'),')\n',
           'Subtes: ',r$subtest,' | Anchor: ',length(r$anchors),'\n',
           'Proporsi sampel: ',r$sample_pct,'% | Model: Rasch (TAM)\n',
           if(!is.null(r$concurrent$error)) paste0('Concurrent: ',r$concurrent$error,'\n') else '',
           if(!is.null(r$fixed$error)) paste0('Fixed: ',r$fixed$error,'\n') else '',
           if(!is.null(r$ctt$error)) paste0('CTT: ',r$ctt$error) else '',
           '\nMissing X: ',round(100*mean(is.na(r$X)),2),'%; Y: ',round(100*mean(is.na(r$Y)),2),'%')
  })
  safe_result <- reactive({
    req(input$run>0)
    tryCatch(result(),error=function(e){validate(need(FALSE,conditionMessage(e)))})
  })
  selected_params <- reactive({
    r <- safe_result()
    if(input$approach=='Concurrent') {
      validate(need(!is.null(r$concurrent$Y),paste('Concurrent gagal:',r$concurrent$error)))
      # Geser hasil concurrent ke baseline hasil kalibrasi terpisah X.
      # Hasil X ditampilkan tetap, sedangkan Y digeser sebesar selisih
      # rata-rata kesukaran anchor (linking lokasi Rasch).
      b_align <- mean(r$parX$b[match(r$anchors,r$parX$item)] -
                      r$concurrent$X$b[match(r$anchors,r$concurrent$X$item)])
      py <- r$concurrent$Y
      py$b <- py$b + b_align
      list(x=r$parX,y=py)
    } else if(input$approach=='Fixed') {
      validate(need(!is.null(r$fixed$par),paste('Fixed gagal:',r$fixed$error)))
      list(x=r$parX,y=r$fixed$par)
    } else if(input$approach=='Linear') list(x=r$parX,y=r$transformed[[input$method]])
    else list(x=r$parX,y=r$parY)
  })
  # Hanya X sebagai baseline; Y sebelum dan Y setelah equating.
  param_data <- reactive({
    r <- safe_result()
    items <- sort(unique(c(r$parX$item,r$parY$item)))
    bx <- r$parX$b[match(items,r$parX$item)]
    by <- r$parY$b[match(items,r$parY$item)]
    yp <- NULL
    if (input$approach == 'Concurrent' && !is.null(r$concurrent$error)) {
      yp <- NULL
    } else if (input$approach == 'Fixed' && !is.null(r$fixed$error)) {
      yp <- NULL
    } else if (input$approach != 'CTT') {
      yp <- selected_params()$y
    }
    by_after <- if(is.null(yp)) rep(NA_real_,length(items)) else yp$b[match(items,yp$item)]
    data.frame(
      Butir=items,
      Jenis=ifelse(items %in% r$anchors,'Anchor','Non-anchor (posisi sama)'),
      X_baseline=bx,
      Y_sebelum=by,
      Y_setelah=by_after,
      Delta_sebelum=by-bx,
      Delta_setelah=by_after-bx,
      stringsAsFactors=FALSE,check.names=FALSE
    )
  })
  output$paramTable <- renderDT({
    d <- param_data()
    if(input$anchor_only) d <- d[d$Jenis=='Anchor',,drop=FALSE]
    show_table(d,pageLength=20)
  })
  diff_data <- reactive({
    d <- param_data()
    if(input$anchor_only) d <- d[d$Jenis=='Anchor',,drop=FALSE]
    d[is.finite(d$X_baseline) & is.finite(d$Y_sebelum),,drop=FALSE]
  })
  output$paramPlot <- renderPlot({
    d <- diff_data()
    validate(need(nrow(d)>0,'Tidak ada parameter yang dapat dibandingkan.'))
    dp <- rbind(
      data.frame(X=d$X_baseline,Y=d$Y_sebelum,Jenis=d$Jenis,Tahap='Sebelum equating'),
      data.frame(X=d$X_baseline,Y=d$Y_setelah,Jenis=d$Jenis,Tahap='Setelah equating')
    )
    dp <- dp[is.finite(dp$X)&is.finite(dp$Y),,drop=FALSE]
    validate(need(nrow(dp)>0,'Parameter hasil equating tidak tersedia.'))
    ggplot(dp,aes(X,Y,color=Jenis))+
      geom_point(size=2,alpha=.75)+
      geom_abline(slope=1,intercept=0,linetype='dashed',colour='grey35')+
      facet_wrap(~Tahap)+theme_minimal()+
      labs(x='Kesukaran X (baseline)',y='Kesukaran Y',
           subtitle='Perbandingan butir identik secara substantif hanya berlaku untuk anchor')
  })
  output$diffTable <- renderDT({show_table(diff_data(),pageLength=20)})

  # Statistik deskriptif parameter butir; tidak menyamakan non-anchor.
  stat_pair <- function(x, y, ids) {
    ix <- match(ids,x$item); iy <- match(ids,y$item)
    valid <- !is.na(ix)&!is.na(iy)
    bx <- x$b[ix[valid]]; by <- y$b[iy[valid]]
    c(Mx=mean(bx),My=mean(by),SDx=sd(bx),SDy=sd(by),
      delta=mean(by)-mean(bx),
      r_par=if(length(bx)>2L && sd(bx)>0 && sd(by)>0) cor(bx,by) else NA_real_)
  }
  # Statistik anchor sebelum equating.
  output$statAnchor <- renderDT({
    r <- safe_result()
    v <- stat_pair(r$parX,r$parY,r$anchors)
    d <- data.frame(group=paste(r$formX,'vs',r$formY),
                    Mx=v['Mx'],My=v['My'],SDx=v['SDx'],SDy=v['SDy'],
                    delta=v['delta'],r_par=v['r_par'],row.names=NULL)
    show_table(d,dom='t')
  })
  # Semua butir: M dan SD untuk masing-masing form, tanpa pemadanan non-anchor.
  output$statAll <- renderDT({
    r <- safe_result()
    d <- data.frame(group=paste(r$formX,'vs',r$formY),
                    Mx=mean(r$parX$b),My=mean(r$parY$b),
                    SDx=sd(r$parX$b),SDy=sd(r$parY$b))
    show_table(d,dom='t')
  })
  output$statAfter <- renderDT({
    r <- safe_result()
    validate(need(input$approach!='CTT','CTT menyetarakan skor total; tidak menghasilkan parameter butir setelah equating.'))
    py <- selected_params()$y
    before <- stat_pair(r$parX,r$parY,r$anchors)
    after <- stat_pair(r$parX,py,r$anchors)
    d <- data.frame(group=paste(r$formX,'vs',r$formY),
                    Mx=after['Mx'],My_sebelum=before['My'],My_setelah=after['My'],
                    SDx=after['SDx'],SDy=after['SDy'],
                    delta=after['delta'],r_par=after['r_par'],row.names=NULL)
    show_table(d,dom='t')
  })
  # Skor total mentah, bukan parameter IRT.
  output$rawPlot <- renderPlot({
    r <- safe_result()
    d <- rbind(data.frame(score=score_display(r$X),form=r$formX),
               data.frame(score=score_display(r$Y),form=r$formY))
    d <- d[is.finite(d$score),,drop=FALSE]
    validate(need(nrow(d)>0,'Tidak ada skor respons yang dapat ditampilkan.'))
    ggplot(d,aes(score,fill=form))+
      geom_histogram(binwidth=1,boundary=-.5,colour='white')+
      facet_wrap(~form,ncol=1,scales='free_y')+
      scale_x_continuous(breaks=pretty(d$score,n=12))+
      theme_minimal()+labs(x='Skor total mentah',y='Jumlah peserta',fill='Form')
  })
  output$rawTable <- renderDT({
    r <- safe_result()
    sx <- score_display(r$X); sy <- score_display(r$Y)
    d <- data.frame(form=c(r$formX,r$formY),N=c(nrow(r$X),nrow(r$Y)),
                    Butir=c(ncol(r$X),ncol(r$Y)),
                    Mean=c(mean(sx,na.rm=TRUE),mean(sy,na.rm=TRUE)),
                    SD=c(sd(sx,na.rm=TRUE),sd(sy,na.rm=TRUE)),
                    Min=c(min(sx,na.rm=TRUE),min(sy,na.rm=TRUE)),
                    Max=c(max(sx,na.rm=TRUE),max(sy,na.rm=TRUE)))
    show_table(d,dom='t')
  })

  icc_items_selected <- reactive({
    r <- safe_result()
    if (input$icc_items == 'anchor') {
      intersect(r$anchors, intersect(r$parX$item, r$parY$item))
    } else {
      intersect(r$parX$item, r$parY$item)
    }
  })
  icc_long <- function(parX, parY, items, labelX, labelY, theta=seq(-4,4,length.out=161)) {
    items <- intersect(items, intersect(parX$item, parY$item))
    validate(need(length(items)>0,'Tidak ada item yang dapat diplot.'))
    make_one <- function(par, form_label) {
      pp <- par[match(items, par$item), , drop=FALSE]
      pr <- prob(theta, pp)
      data.frame(theta = rep(theta, times=length(items)),
                 p = as.vector(pr),
                 item = rep(items, each=length(theta)),
                 form = form_label)
    }
    rbind(make_one(parX, labelX), make_one(parY, labelY))
  }
  plot_icc_faceted <- function(dat, ttl, subttl) {
    ggplot(dat, aes(theta, p, color=form)) +
      geom_line(linewidth=.85) +
      facet_wrap(~item, ncol=5) +
      theme_bw() +
      theme(
        strip.background = element_rect(fill='grey90', colour='grey50'),
        panel.grid.minor = element_line(colour='grey92'),
        panel.grid.major = element_line(colour='grey88')
      ) +
      labs(title=ttl, subtitle=subttl, x='Theta', y='P(theta)', color='Kelompok')
  }
  output$iccPrePlot <- renderPlot({
    r <- safe_result()
    items <- icc_items_selected()
    d <- icc_long(r$parX, r$parY, items, r$formX, r$formY)
    plot_icc_faceted(
      d,
      ttl = paste0('ICC sebelum equating: ', r$formX, ' vs ', r$formY),
      subttl = if (input$icc_items == 'anchor') 'Hanya item anchor' else 'Semua item yang sama pada kedua form'
    )
  })
  output$iccPostPlot <- renderPlot({
    r <- safe_result()
    validate(need(input$approach != 'CTT','Tab ini hanya untuk pendekatan IRT (Concurrent, Fixed, atau Transformasi Linier).'))
    p <- selected_params()
    items <- icc_items_selected()
    d <- icc_long(p$x, p$y, items, r$formX, r$formY)
    ttl_method <- if (input$approach == 'Linear') paste0('Linear - ', input$method) else input$approach
    plot_icc_faceted(
      d,
      ttl = paste0('ICC setelah equating: ', r$formX, ' vs ', r$formY),
      subttl = paste('Pendekatan:', ttl_method)
    )
  })
  # Konversi IRT expected-score: inversi TCC X untuk TCC Y
  conversion <- reactive({
    r <- safe_result()
    if(input$approach=='CTT') {
      validate(need(!is.null(r$ctt$tab),paste('CTT gagal:',r$ctt$error)))
      return(r$ctt$tab)
    }
    p <- selected_params()
    th <- seq(-8,8,length.out=2001)
    tX <- rowSums(prob(th,p$x)); tY <- rowSums(prob(th,p$y))
    # peta skor harapan Y ke skor harapan X pada theta sama
    target <- 0:ncol(r$Y)
    theta_at_score <- approx(x=tY,y=th,xout=pmin(max(tY),pmax(min(tY),target)),ties='ordered',rule=2)$y
    data.frame(score_Y=target,score_X_equated=approx(th,tX,xout=theta_at_score,rule=2)$y)
  })
  output$conversionTable <- renderDT({show_table(conversion(),pageLength=12)})
  output$equatedPlot <- renderPlot({
    r <- safe_result(); cv <- conversion()
    scoreY <- score_display(r$Y)
    mapped <- approx(cv$score_Y,cv$score_X_equated,xout=scoreY,rule=2)$y
    d <- rbind(data.frame(score=score_display(r$X),series='X: skor mentah'),
               data.frame(score=scoreY,series='Y: sebelum'),
               data.frame(score=mapped,series='Y: setelah'))
    ggplot(d,aes(score,fill=series))+geom_histogram(bins=30,alpha=.7)+
      facet_wrap(~series,scales='free_y',ncol=1)+theme_minimal()+
      labs(x='Skor total (skala X setelah equating)',y='Frekuensi',fill=NULL)
  })
  tcc_figure <- function(px,py,labelX,labelY,title,subtitle) {
    th <- seq(-4,4,length.out=201)
    d <- rbind(data.frame(theta=th,score=rowSums(prob(th,px)),form=labelX),
               data.frame(theta=th,score=rowSums(prob(th,py)),form=labelY))
    ggplot(d,aes(theta,score,color=form))+
      geom_line(linewidth=1)+theme_minimal()+
      labs(title=title,subtitle=subtitle,x='Theta',
           y='Expected total score (TCC)',color='Form')
  }
  output$tccBefore <- renderPlot({
    r <- safe_result()
    tcc_figure(r$parX,r$parY,r$formX,r$formY,
               'TCC sebelum equating',
               'Kalibrasi terpisah; sumbu theta belum memiliki metrik bersama')
  })
  output$tccAfter <- renderPlot({
    r <- safe_result()
    validate(need(input$approach!='CTT',
      'CTT tidak mengubah parameter IRT; plot TCC setelah equating tidak tersedia.'))
    p <- selected_params()
    lbl <- if(input$approach=='Linear') paste('Linear',input$method) else input$approach
    tcc_figure(p$x,p$y,r$formX,r$formY,
               'TCC setelah equating',paste('Metode:',lbl,'— pada metrik X'))
  })
  output$linkTable <- renderDT({show_table(safe_result()$links,dom='t')})
  output$linkPlot <- renderPlot({
    r <- safe_result(); d <- r$links
    d2 <- rbind(data.frame(method=d$method,parameter='A',value=d$A),
                data.frame(method=d$method,parameter='B',value=d$B))
    ggplot(d2,aes(method,value,fill=method))+geom_col()+facet_wrap(~parameter,scales='free_y')+
      theme_minimal()+labs(x='Metode',y='Konstanta linking')+guides(fill='none')
  })
  output$summaryTable <- renderDT({
    r <- safe_result()
    d <- data.frame(approach=c('Concurrent','Fixed','Linear MM','Linear MS','Linear SL','Linear HB','CTT'),
      status=c(if(is.null(r$concurrent$error)) 'Berhasil' else r$concurrent$error,
               if(is.null(r$fixed$error)) 'Berhasil' else r$fixed$error,
               rep('Berhasil',4),if(is.null(r$ctt$error)) 'Berhasil' else r$ctt$error))
    show_table(d,dom='t')
  })
  output$missingTable <- renderDT({
    r <- safe_result()
    z <- lapply(list(r$X,r$Y),function(d) c(N_awal=attr(d,'original_n'),baris_NA_dibuang=attr(d,'removed_all_na'),N_valid=attr(d,'valid_n'),persen_sampel=attr(d,'sample_pct'),N_dipakai=nrow(d),butir=ncol(d),missing=sum(is.na(d)),persen_missing=100*mean(is.na(d)),peserta_lengkap=sum(complete.cases(d))))
    show_table(data.frame(form=c(input$formX,input$formY),do.call(rbind,z),row.names=NULL),dom='t')
  })
  output$notes <- renderText({
    r <- safe_result()
    paste('Subtes:',r$subtest,'\nAnchor:',paste(r$anchors,collapse=', '),
          '\nPerhatian: untuk item non-anchor, posisi kolom yang sama tidak membuktikan identitas item yang sama.',
          '\nPada model Rasch A=1. CTT menghasilkan konversi skor total, bukan parameter butir.',
          '\nPilihan form dan persentase sampel baru berlaku setelah tombol Jalankan equating diklik.')
  })
}
shinyApp(ui,server)
