# EQUATING BERPASANGAN: form201 (X) vs form202..210 (Y)
# Instal sekali: install.packages(c('mirt','equateIRT'))
library(mirt)

folder <- 'D:/OneDrive/Github/Learning via R/dataset equating'
files <- file.path(folder, sprintf('form%d.csv', 201:205))
stopifnot(all(file.exists(files)))
model <- 'Rasch'  # ganti menjadi '2PL' bila diperlukan
anchor <- sprintf('%02d', c(2,5,9,13,16,20,24,29,33,38,43,47,48,50,52,66,67,68,73,83,84,85,86))
outdir <- file.path(folder, 'hasil_equating')
dir.create(outdir, showWarnings = FALSE, recursive = TRUE)

read_form <- function(path) {
  z <- read.csv(path, header = FALSE, check.names = FALSE, na.strings = c('', 'NA'))
  if (ncol(z) < max(as.integer(anchor))) stop('Kolom kurang dari nomor anchor tertinggi: ', path)
  names(z) <- sprintf('%02d', seq_len(ncol(z)))
  z[] <- lapply(z, function(x) suppressWarnings(as.numeric(as.character(x))))
  v <- unlist(z, use.names = FALSE)
  if (any(!is.na(v) & !v %in% c(0,1))) stop('Data bukan respons 0/1: ', path)
  z
}

# Ambil a dan b (IRT) dari model satu kelompok
get_ab <- function(fit) {
  p <- as.data.frame(mirt::coef(fit, IRTpars = TRUE, simplify = TRUE)$items)
  data.frame(item = rownames(p), a = if ('a' %in% names(p)) p$a else 1,
             b = p$b, row.names = NULL)
}

# Kembalikan tabel parameter untuk multipleGroup, grup 1=X dan 2=Y
get_mg_ab <- function(fit, group) {
  p <- as.data.frame(mirt::coef(fit, IRTpars = TRUE, simplify = TRUE)[[group]]$items)
  data.frame(item = rownames(p), a = if ('a' %in% names(p)) p$a else 1,
             b = p$b, row.names = NULL)
}

# Metode transformasi Y -> X: bX = A*bY+B; aX=aY/A
# MM: rasio rata-rata diskriminasi; MS: rasio simpangan baku b
# SL/HB: minimisasi fungsi respons pada grid theta skala X.
link_ab <- function(px, py, method, grid = seq(-4,4,length.out=121)) {
  ix <- match(anchor, px$item); iy <- match(anchor, py$item)
  x <- px[ix, ]; y <- py[iy, ]
  if (anyNA(ix) || anyNA(iy)) stop('Anchor tidak lengkap pada hasil kalibrasi')
  if (method == 'MM') {
    A <- if (model == 'Rasch') 1 else mean(y$a)/mean(x$a)
    B <- mean(x$b) - A*mean(y$b)
  } else if (method == 'MS') {
    A <- if (model == 'Rasch') 1 else sd(x$b)/sd(y$b)
    B <- mean(x$b) - A*mean(y$b)
  } else {
    pfun <- function(theta, a, b) plogis(outer(theta, b, '-') *
                                         matrix(a, length(theta), length(a), byrow=TRUE))
    target <- pfun(grid, x$a, x$b)
    objective <- function(z) {
      A <- if (model == 'Rasch') 1 else exp(z[1])
      B <- if (model == 'Rasch') z[1] else z[2]
      candidate <- pfun(grid, y$a/A, A*y$b+B)
      if (method == 'SL') sum((rowSums(target)-rowSums(candidate))^2)
      else sum((target-candidate)^2)
    }
    start <- if (model == 'Rasch') 0 else c(0, mean(x$b)-mean(y$b))
    opt <- optim(start, objective, method='BFGS')
    A <- if (model == 'Rasch') 1 else exp(opt$par[1])
    B <- if (model == 'Rasch') opt$par[1] else opt$par[2]
  }
  data.frame(method=method, A=A, B=B)
}

# Fixed item calibration: kunci parameter internal (a1,d) anchor X pada Y
fixed_calibrate <- function(Xfit, Y) {
  pX <- as.data.frame(mirt::coef(Xfit, simplify=TRUE)$items)
  vals <- mirt::mirt(Y, 1, itemtype=model, pars='values', verbose=FALSE)
  for (id in anchor) {
    for (par in c('a1','d')) {
      row <- which(vals$item == id & vals$name == par)
      if (length(row) != 1 || !par %in% names(pX)) stop('Parameter fixed tidak ditemukan: ', id, '/', par)
      vals$value[row] <- pX[id,par]
      vals$est[row] <- FALSE
    }
  }
  mirt::mirt(Y, 1, itemtype=model, pars=vals, verbose=FALSE,
             technical=list(NCYCLES=1000))
}

# CTT: chained linear equating Y total -> anchor Y -> anchor X -> X total
ctt_equate <- function(X,Y) {
  tx <- rowSums(X,na.rm=FALSE); ty <- rowSums(Y,na.rm=FALSE)
  ax <- rowSums(X[,anchor,drop=FALSE],na.rm=FALSE)
  ay <- rowSums(Y[,anchor,drop=FALSE],na.rm=FALSE)
  # Mean/SD anchor dipakai sebagai penghubung; kelompok tidak diasumsikan ekuivalen
  slope <- sd(tx,na.rm=TRUE)/sd(ax,na.rm=TRUE) * sd(ay,na.rm=TRUE)/sd(ty,na.rm=TRUE)
  intercept <- mean(tx,na.rm=TRUE) - slope*mean(ty,na.rm=TRUE) +
    sd(tx,na.rm=TRUE)/sd(ax,na.rm=TRUE) * (mean(ax,na.rm=TRUE)-mean(ay,na.rm=TRUE))
  # CTT linear chaining: X-score = mean(Tx) + (sd(Tx)/sd(Ax))*
  # [mean(Ay)-mean(Ax) + (sd(Ay)/sd(Ty))*(Ty-mean(Ty))]
  # Tanda intercept disesuaikan dengan ekspresi di atas:
  intercept <- mean(tx,na.rm=TRUE) + sd(tx,na.rm=TRUE)/sd(ax,na.rm=TRUE)*
    (mean(ay,na.rm=TRUE)-mean(ax,na.rm=TRUE)) - slope*mean(ty,na.rm=TRUE)
  data.frame(score_Y=0:ncol(Y), score_X_equated=intercept+slope*(0:ncol(Y)))
}

X <- read_form(files[1])
if (!all(anchor %in% names(X))) stop('Anchor tidak ada pada form201')
fitX <- mirt::mirt(X, 1, itemtype=model, verbose=FALSE)
parX <- get_ab(fitX)
all_summary <- list()
all_parameters <- list()
all_ctt <- list()
all_errors <- list()

for (i in 2:length(files)) {
  formY <- tools::file_path_sans_ext(basename(files[i]))
  message('\n=== form201 vs ', formY, ' ===')
  tryCatch({
    Y <- read_form(files[i])
    if (!all(anchor %in% names(Y))) stop('Anchor tidak lengkap')
    # Gabungkan kolom unik: anchor tetap ID yang sama;
    # item nonanchor diberi nama berbeda agar tidak dipaksakan sama.
    nonX <- setdiff(names(X),anchor); nonY <- setdiff(names(Y),anchor)
    namesX <- c(anchor, paste0('X_',nonX))
    namesY <- c(anchor, paste0('Y_',nonY))
    Xc <- X[,c(anchor,nonX),drop=FALSE]; names(Xc) <- namesX
    Yc <- Y[,c(anchor,nonY),drop=FALSE]; names(Yc) <- namesY
    common <- union(names(Xc),names(Yc))
    pad <- function(d) { z <- matrix(NA_real_,nrow(d),length(common),dimnames=list(NULL,common)); z[,names(d)] <- as.matrix(d); as.data.frame(z) }
    joined <- rbind(pad(Xc),pad(Yc))
    grp <- factor(c(rep('X',nrow(X)),rep('Y',nrow(Y))),levels=c('X','Y'))
    # 1. Concurrent: anchor invariant, distribusi theta grup Y bebas
    fitConcurrent <- mirt::multipleGroup(joined, 1, group=grp, itemtype=model,
      invariance=c(anchor,'free_means','free_var'), verbose=FALSE)
    pcX <- get_mg_ab(fitConcurrent,'X'); pcY <- get_mg_ab(fitConcurrent,'Y')
    pcX$item <- sub('^X_','',pcX$item); pcY$item <- sub('^Y_','',pcY$item)
    pcX <- pcX[!grepl('^Y_',pcX$item),]; pcY <- pcY[!grepl('^X_',pcY$item),]
    pcX <- pcX[pcX$item %in% names(X),]; pcY <- pcY[pcY$item %in% names(Y),]
    pcX$approach <- 'Concurrent X'; pcY$approach <- 'Concurrent Y'
    # 2. Fixed anchor X -> Y
    fitFixed <- fixed_calibrate(fitX,Y)
    pf <- get_ab(fitFixed); pf$approach <- 'Fixed Y'
    # 3. Separate calibration + MM/MS/SL/HB
    fitY <- mirt::mirt(Y,1,itemtype=model,verbose=FALSE)
    parY <- get_ab(fitY)
    links <- do.call(rbind,lapply(c('MM','MS','SL','HB'), function(m) link_ab(parX,parY,m)))
    transformed <- do.call(rbind,lapply(seq_len(nrow(links)),function(j) {
      z <- parY; z$a <- z$a/links$A[j]; z$b <- links$A[j]*z$b+links$B[j]
      z$approach <- paste0('Linear ',links$method[j]); z
    }))
    # 4. CTT score conversion
    eqCTT <- ctt_equate(X,Y)
    eqCTT$form <- formY
    # Simpan tabel
    params <- rbind(transform(parX,approach='Baseline X'),pcX,pcY,pf,transformed)
    params$form <- formY
    links$form <- formY
    links <- links[,c('form','method','A','B')]
    all_summary[[formY]] <- links
    all_parameters[[formY]] <- params
    all_ctt[[formY]] <- eqCTT
    write.csv(links,file.path(outdir,paste0(formY,'_linking.csv')),row.names=FALSE)
    write.csv(params,file.path(outdir,paste0(formY,'_parameters.csv')),row.names=FALSE)
    write.csv(eqCTT,file.path(outdir,paste0(formY,'_CTT.csv')),row.names=FALSE)
    print(links,row.names=FALSE)
  }, error=function(e) {
    message('GAGAL ',formY,': ',conditionMessage(e))
    all_errors[[formY]] <<- conditionMessage(e)
  })
}
if(length(all_summary)) write.csv(do.call(rbind,all_summary),file.path(outdir,'RINGKASAN_LINKING.csv'),row.names=FALSE)
if(length(all_parameters)) write.csv(do.call(rbind,all_parameters),file.path(outdir,'SEMUA_PARAMETER.csv'),row.names=FALSE)
if(length(all_ctt)) write.csv(do.call(rbind,all_ctt),file.path(outdir,'SEMUA_CTT.csv'),row.names=FALSE)
if(length(all_errors)) write.csv(data.frame(form=names(all_errors),error=unlist(all_errors)),file.path(outdir,'ERROR.csv'),row.names=FALSE)
message('Selesai. Output: ',outdir)
