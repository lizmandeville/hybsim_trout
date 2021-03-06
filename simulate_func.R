simulate.hyb <- function(nind.start = 100, prop.sp1=0.5, n.generation=10, makeplot=TRUE, printoutput=TRUE){

## Set up matrix for keeping track of q, or proportion of ancestry
## Since q1+q2=1, will only keep track of q for 1 category
    q.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind.start)
    q.allgen[1,] <- c(rep(1,round(prop.sp1*nind.start, digits=0)), rep(0,(nind.start-round(prop.sp1*nind.start, digits=0))))

    qfromQ.allgen <-matrix(data=NA, nrow=(n.generation+1), ncol=nind.start)
    qfromQ.allgen[1,] <- q.allgen[1,]

    Q11.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind.start)
    Q11.allgen[1,] <- c(rep(1,round(prop.sp1*nind.start, digits=0)), rep(0,(nind.start-round(prop.sp1*nind.start, digits=0))))

    Q12.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind.start)
    Q12.allgen[1,] <- rep(0, nind.start)

    Q21.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind.start)
    Q21.allgen[1,] <- rep(0, nind.start)

    Q22.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind.start)
    Q22.allgen[1,] <- c(rep(0,round(prop.sp1*nind.start, digits=0)), rep(1,(nind.start-round(prop.sp1*nind.start, digits=0))))

    Qinter.allgen <- matrix(data=NA, nrow=(n.generation+1), ncol=nind.start)
    Qinter.allgen[1,] <- rep(0, nind.start)

    if(makeplot == TRUE){
        pdf(paste("qQ_",n.generation,"gen_",prop.sp1,"sp1_",nind.start,"ind.pdf",sep=""))
##qfromQ <- Q11.allgen+((Q12.allgen+Q21.allgen)/2)
        par(mfrow=c(ceiling(n.generation/4),4))
    }

    for(i in 2:(n.generation+1)){

    ## sample a vector of integers 1:nind, then use that vector to index q.allgen and Q.allgen so we can get matching q and Q
    ## sampling WITH REPLACEMENT, i.e. the fact that an individual has already produced one offspring does not make it ineligible to produce another
        randvec1 <- sample(1:nind.start, nind.start, replace=T)
        randvec2 <- sample(1:nind.start, nind.start, replace=T)

        ## sample from prev. generation for new offspring, assume random mating
        ## use random vectors to index both q and Q
        ## once I'm convinced that I'm getting q from Q correctly, will drop this part
        ##parent1.q <- q.allgen[i-1,randvec1]
        ##parent2.q <- q.allgen[i-1,randvec2]

        ## get q for offspring - this will go away too
        ##q.allgen[i,] <- (parent1.q+parent2.q)/2

        parent1.Q11 <- Q11.allgen[i-1,randvec1]
        parent1.Q12 <- Q12.allgen[i-1,randvec1]
        parent1.Q21 <- Q21.allgen[i-1,randvec1]
        parent1.Q22 <- Q22.allgen[i-1,randvec1]

        parent2.Q11 <- Q11.allgen[i-1,randvec2]
        parent2.Q12 <- Q12.allgen[i-1,randvec2]
        parent2.Q21 <- Q21.allgen[i-1,randvec2]
        parent2.Q22 <- Q22.allgen[i-1,randvec2]

        Q11.allgen[i,] <- (parent1.Q11+parent1.Q12)*(parent2.Q11+parent2.Q12)
        Q12.allgen[i,] <- (((parent1.Q11+parent1.Q12)*(parent2.Q21+parent2.Q22))+((parent1.Q21+parent1.Q22)*(parent2.Q11+parent2.Q12)))/2 ## Need to do this to retain symmetry of Q12,Q21
        Q21.allgen[i,] <- Q12.allgen[i,]
        Q22.allgen[i,] <- (parent1.Q21+parent1.Q22)*(parent2.Q21+parent2.Q22)

        qfromQ.allgen[i,] <- Q11.allgen[i,]+((Q12.allgen[i,]+Q21.allgen[i,])/2)

        Qinter.allgen[i,] <- Q12.allgen[i,]+Q21.allgen[i,]

        if(makeplot == TRUE){
            plot(qfromQ.allgen[i,], Qinter.allgen[i,], main=paste("Gen.", i-1, sep=" "), xlab="q", ylab="Q", xlim=c(0,1), ylim=c(0,1), type="n", axes=F)
    ##plot(q.allgen[i,], qfromQ.allgen[i,], main=paste("Gen.", i-1, sep=" "), xlab="q", ylab="qfromQ", xlim=c(0,1), ylim=c(0,1))
            axis(1, at=c(0,0.5,1), labels=c(0,0.5,1))
            axis(2, at=c(0,0.5,1), labels=c(0,0.5,1))
            arrows(0,0,0.5,1, length=0, col="gray90")
            arrows(0.5,1,1,0, length=0, col="gray90")
            points(qfromQ.allgen[i,], Qinter.allgen[i,], col="gray45")
        }
    }
    if(makeplot == TRUE){
        dev.off()

    }

    if(printoutput == TRUE){
        return(list(qfromQ.allgen, Qinter.allgen))
    }
}

## This is not realistic because in the real NF, they kept adding individuals from parental species by stocking. Simulate this?
