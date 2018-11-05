
n.draws <- 100000
n.thin <- 10
n.burn <- 2000
datadir <- "data"
outdir <- "out"

library(magrittr)
library(dplyr)

library(lubridate)
library(readr)
library(nimble)
library(coda)
library(reshape2)
library(scales)

source("color.R")

all.fn <- list.files(datadir, full.names=TRUE)
slugs <- all.fn %>% basename %>% gsub("\\.csv", "", .)

n <- length(all.fn)

should.plot.trace <- FALSE

all.out <- list()
for (ix in 1:length(all.fn)) {
    fn.dat <- all.fn[ix]
    slug <- slugs[ix]
    print.(sprintf("iteration %d of %d - race %s", ix, n, slug), newline="") 

    dat.raw <- read_csv(fn.dat)
    colnames(dat.raw) <-
        colnames(dat.raw) %>% 
        gsub("[()]", "", .) %>%
        gsub(" ", ".", .)

    dat <- dat.raw %>%
        filter(Sample != "LV") %>% 
        mutate(N = as.numeric(gsub(" ([RL]V|A)", "", Sample)),
               VoterType = gsub("^.* ", "", Sample),
               Start = gsub(" .*$", "", Date),
               StartDate = mdy(sprintf("%s/%d", Start, Year)),
               End = gsub("^.* ", "", Date),
               EndDate = mdy(sprintf("%s/%d", End, Year))) %>%
        select(-Spread, -Sample) %>%
        arrange(desc(EndDate), desc(StartDate))
    dat$i <- 1:nrow(dat)

    other <- dat %>% select(-Date, -Poll,
                            -N, -VoterType,
                            -Start, -End, -StartDate, -Year)
    if ("MoE" %in% colnames(other)) {
        other <- other %>% select(-MoE)
    }
    other$Other <- 100 - rowSums(select(other, -EndDate, -i))
    no.date <- !(colnames(other) %in% c("EndDate", "i"))
    other[no.date] <- other[no.date]/100

    series.names <- colnames(select(other, -EndDate, -i))
    dat <- dat[!colnames(dat) %in% series.names] %>% left_join(other)

    date.range <- range(dat$EndDate)
    date.seq <- seq(date.range[1], date.range[2],
                    by="1 d")
    T <- length(date.seq)

    dates <- data.frame(Date=date.seq,
                        Index=1:T,
                        stringsAsFactors=FALSE)
    polls <- dat %>% left_join(dates, by=c("EndDate" = "Date"))

    x <- polls[,c("Index", "N", series.names)]
    for (sn in series.names) {
        x[sprintf("%s.N", sn)] <- x[,sn]*x$N
    }
    series.names.N <- sprintf("%s.N", series.names)
    x <- x[,!colnames(x) %in% series.names]
    x <- x %>%
        group_by(Index) %>%
        summarize_at(vars(-N, -Index), sum)
    x$N <- select(x, -Index) %>% rowSums
    colnames(x) <- gsub("\\.N$", "", colnames(x))
    x[series.names] <- x[series.names] / x$N

    r <- length(series.names)-1
    N <- nrow(x)
    ny <- x %>% group_by(Index) %>% summarize(ny=n())
    
    all.ny <- dates %>% left_join(ny, by="Index") %>% mutate(ny = replace(ny, is.na(ny), 0))
    ny.0 <- all.ny %>% filter(ny == 0)
    ny.0$ix0 <- 1:nrow(ny.0)
    ny0 <- nrow(ny.0)
    ny.0 <- all.ny %>% left_join(ny.0, by=c("Index", "Date"))

    all.x0 <- x %>% right_join(dates)
    all.x <- all.x0 %>% 
        group_by(Index, Date) %>%
        summarize_all(funs(mean)) %>%
        mutate(N = round(N)) %>% 
        ungroup %>%
        left_join(all.x0)

    y <- sapply(1:nrow(x), function(i) {
        row <- x[i,]
        o <- unlist(row[series.names] * row$N)
        o <- round(o[names(o) != "Other"])
        o <- c(o, Other = row$N - sum(round(o)))
        o
    }) %>% t

    all.y <- sapply(1:nrow(all.x), function(i) {
        row <- all.x[i,]
        o <- unlist(row[series.names] * row$N)
        o <- round(o[names(o) != "Other"])
        o <- c(o, Other = row$N - sum(round(o)))
        o
    }) %>% t

    p_to_mvn <- function(mat) {
        apply(mat, 1,
              function(x) {
                  p. <- tail(x, 1)
                  log(head(x, -1) / p.)
              }) %>% t
    }
    eta <- p_to_mvn(all.y) %>% t
    all.y <- t(all.y)
    all.y[,all.ny$ny > 1] <- NA

    inputs <- list(
        const=list(r = r,
                   N = N,
                   T = T,
                   b0 = 0.001,
                   a0 = 0.001,
                   C0 = 1,
                   k = x$Index,
                   ss = x$N,
                   all.y = all.y,
                   ny = all.ny$ny),
        data=list(filtered = array(0, c(r, T)),
                  m0 = rep(0, r),
                  B0 = 1,
                  B = rep(1, T),
                  C = rep(1, T),
                  R = rep(1, T),
                  h = array(0, c(r,T)),
                  H = array(1, c(r,r,T)),
                  h0 = rep(0, r),
                  H0 = array(1, c(r,r)),
                  y = y,
                  pi = array(1/(r+1), c(r+1,T)),
                  eta = eta),
        inits=list(beta=array(0, c(r, T)),
                   v = 1,
                   w = 1,
                   beta0 = rep(0, r))
    )

    get.filtered <- nimbleFunction(
        run = function(ny = integer(0), R = double(0), v = double(0), m = double(1), eta = double(1), r = integer(0)) {
            returnType(double(1))
            if (ny == 0) {
                return(m)
            }
            A <- R / (R+v)
            filtered <- m[1:r] + A*(eta[1:r]-m[1:r])
            return(filtered)
        }
    )
    get.C <- nimbleFunction(
        run = function(R = double(0), v = double(0)) {
            returnType(double(0))
            A <- R / (R+v)
            C <- R - A*A*(R+v)
            return(C)
        }
    )
    p.to.mvn <- nimbleFunction(
        run = function(x = double(1)) {
            returnType(double(1))
            p. <- x[length(x)]
            x. <- x[-length(x)]
            return(log(x. / p.))
        })
    mvn.to.p <- nimbleFunction(
        run = function(x = double(1)) {
            returnType(double(1))
            vv <- c(exp(x), 1)
            return(vv / sum(vv))
        })

    code <- nimbleCode({
        ## Forward-filter
        R[1] <- C0 + w
        C[1] <- get.C(R[1], v)
        filtered[1:r,1] <- get.filtered(ny[1], R[1], v, m0[1:r], eta[1:r,1], r)
        B0 <- C0 / R[1]
        ##
        for (i in 2:T) {
            R[i] <- C[i-1] + w
            filtered[1:r,i] <- get.filtered(ny[i], R[i], v,
                                            filtered[1:r,i-1],
                                            eta[1:r,i],
                                            r)
            B[i-1] <- C[i-1] / R[i]
        }
        ## Backwards sampling
        H[1:r,1:r,T] <- diag(rep(C[T], r))
        beta[1:r,T] ~ dmnorm(mean=filtered[1:r,T],
                             cov=H[1:r,1:r,T])
        pi[1:(r+1),T] <- mvn.to.p(beta[1:r,T])
        for (i in 1:(T-1)) {
            h[1:r,T-i] <- filtered[1:r,T-i] + B[T-i]*(beta[1:r,T-i+1]-filtered[1:r,T-i+1])
            H[1:r,1:r,T-i] <- diag(rep(C[T-i] - B[T-i]*B[T-i]*R[T-i+1], r))
            beta[1:r,T-i] ~ dmnorm(mean=h[1:r,T-i], cov=H[1:r,1:r,T-i])
            ##
            pi[1:(r+1),T-i] <- mvn.to.p(beta[1:r,T-i])
        }
        h0[1:r] <- m0[1:r] + B0*(beta[1:r,1]-filtered[1:r,1])
        H0[1:r,1:r] <- diag(rep(C0 - B0*B0*R[1], r))
        beta0[1:r] ~ dmnorm(mean=h0[1:r], cov=H0[1:r,1:r])
        ##
        w ~ dinvgamma(shape=a0, scale=b0)
        v ~ dinvgamma(shape=a0, scale=b0)
        for (i in 1:N) {
            ## Likelihood
            y[i, 1:(r+1)] ~ dmulti(size=ss[i],
                                   prob=pi[1:(r+1),k[i]])
        }
    })
    ##
    Rmodel <- with(inputs, {
        nimbleModel(code = code, constants = const,
                    data = data, inits = inits)
    })

    RmodelSpec <- configureMCMC(Rmodel, print=TRUE,
                                monitors=c("v", "w", "pi"),
                                thin=n.thin, nburnin=n.burn,
                                useConjugacy=FALSE)
    Rmcmc <- buildMCMC(RmodelSpec, useConjugacy=FALSE)
    Cmodel <- compileNimble(Rmodel, showCompilerOutput=FALSE)
    Cmcmc <- compileNimble(Rmcmc, project=Rmodel, showCompilerOutput=FALSE)

    Cmcmc$run(n.draws)
    m <- as.mcmc(as.matrix(Cmcmc$mvSamples))

    ## if (should.plot.trace) {
    ##     ix.pi.wv <- grepl(paste(c("v","w", sprintf(" %d]", x$Index)), collapse="|"), colnames(m))
    ##     plot(m[,ix.pi.wv])
    ## }

    ix.D <- which(grepl("D$", series.names))
    ix.R <- which(grepl("R$", series.names))
    ix.L <- which(grepl("L$", series.names))
    ix.I <- which(grepl("I$", series.names))
    ix.all <- c(ix.D, ix.R, ix.L, ix.I)
    series.names.cols <- series.names[c(ix.all, which(!grepl("[RDLI]$", series.names)))]

    x.dated <- x %>% left_join(dates)

    p.D.win <- mean(apply(m[,sprintf("pi[%d, %d]", ix.all, T)], 1,
                          function(x) which.max(x) == 1))
    p.R.win <- mean(apply(m[,sprintf("pi[%d, %d]", ix.all, T)], 1,
                          function(x) which.max(x) == 2))
    p.L.win <- mean(apply(m[,sprintf("pi[%d, %d]", ix.all, T)], 1,
                          function(x) which.max(x) == 3))
    p.I.win <- mean(apply(m[,sprintf("pi[%d, %d]", ix.all, T)], 1,
                          function(x) which.max(x) == 4))
    
    cols <- cols_all[c("Blue2", "Red2", "Orange2", "GreenB2", "Purple2")]
    cols.light <- cols_all[c("Blue2", "Red2", "Orange2", "GreenB2", "Purple2")]
    names(cols) <- names(cols.light) <- c("D", "R", "L", "I", "Other")
    ##
    ix.cols <- gsub(".*\\.", "", series.names.cols)
    cols <- cols[ix.cols]
    cols.light <- cols.light[ix.cols]
    names(cols) <- names(cols.light) <- series.names.cols
    ##
    post.mean <- colMeans(m)
    post.sd <- apply(m, 2, sd)
    ##
    pi.hat <- post.mean[grepl("pi\\[", names(post.mean))]
    pi.hat <- data.frame(pi = unname(pi.hat),
                         j = gsub("pi\\[|,.*\\]$", "", names(pi.hat)),
                         i = gsub(".*, |\\]", "", names(pi.hat))) %>% 
        acast(i ~ j, value.var="pi")
    pi.hat <- pi.hat[order(as.numeric(rownames(pi.hat))),]
    colnames(pi.hat) <- series.names
    ##
    pi.sd <- post.sd[grepl("pi\\[", names(post.mean))]
    pi.sd <- data.frame(pi = unname(pi.sd),
                        j = gsub("pi\\[|,.*\\]$", "", names(pi.sd)),
                        i = gsub(".*, |\\]", "", names(pi.sd))) %>% 
        acast(i ~ j, value.var="pi")
    pi.sd <- pi.sd[order(as.numeric(rownames(pi.sd))),]
    colnames(pi.sd) <- series.names
    ##
    post.b <- apply(m, 2, quantile, c(0.025, 0.975))
    ##
    pi.ub <- post.b[2,grepl("pi\\[", names(post.b[2,]))]
    pi.ub <- data.frame(pi = unname(pi.ub),
                        j = gsub("pi\\[|,.*\\]$", "", names(pi.ub)),
                        i = gsub(".*, |\\]", "", names(pi.ub))) %>% 
        acast(i ~ j, value.var="pi")
    pi.ub <- pi.ub[order(as.numeric(rownames(pi.ub))),]
    colnames(pi.ub) <- series.names
    ##
    pi.lb <- post.b[1,grepl("pi\\[", names(post.b[2,]))]
    pi.lb <- data.frame(pi = unname(pi.lb),
                        j = gsub("pi\\[|,.*\\]$", "", names(pi.lb)),
                        i = gsub(".*, |\\]", "", names(pi.lb))) %>% 
        acast(i ~ j, value.var="pi")
    pi.lb <- pi.lb[order(as.numeric(rownames(pi.lb))),]
    colnames(pi.lb) <- series.names
    ##
    fn.out <- sprintf("%s/%s.png", outdir, slug)
    png(fn.out, width=1600, height=800)
    plot(range(all.x$Date), range(pi.lb, pi.ub), type="n",
         main=sprintf("Prob.(Democrat wins) = %.2f", p.D.win))
    for (i in 1:(r+1)) {
        sn <- series.names[i]
        if (sn == "Other") next
        polygon(c(all.x$Date, rev(all.x$Date)),
                c(pi.ub[,sn], rev(pi.lb[,sn])),
                col=alpha(cols.light[sn], .1),
                border=NA)
        ##
        ## lines(polls$EndDate, x[[sn]], lty=1, col=cols[sn])
        points(x.dated$Date, x.dated[[sn]], pch=19, col=cols[sn])
        lines(all.x$Date, pi.hat[,sn], lwd=3, col=cols[sn])
        lines(all.x$Date, pi.ub[,sn], lwd=2, col=cols[sn], lty=2)
        lines(all.x$Date, pi.lb[,sn], lwd=2, col=cols[sn], lty=2)
    }
    legend("topleft", legend=series.names, border=cols[series.names], fill=cols.light[series.names], density=rep(NA, 3),
           bty="n")
    dev.off()

    rm(m, Rmodel, RmodelSpec, Rmcmc, Cmodel, Cmcmc)
    gc()

    all.out[[slug]] <- list(slug=slug, all.x=all.x, polls=polls, pi.lb=pi.lb, pi.ub=pi.ub, p.D.win=p.D.win, post.mean=post.mean, pi.hat=pi.hat, post.sd=post.sd)
    cat("Done.\n")
}
