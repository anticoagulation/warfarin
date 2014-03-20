##########
#
# Functions to calculate various measures of INR control, including
# Time in Therapeutic Range (TTR) following the linear interpolation 
# method of Rosendaal.
#
# Project started 2014-03-05 by Al Ozonoff, Boston Children's Hospital
# 2014-03-17 made changes to accept multiple patient id's
#
# Some of this research was funded by VA grant XXX PI Adam Rose
# 
# Licensed under GNU Public License v3
#
##########






# function draw.inr
# 
# Graph a list of INR values with a single longitudinal plot.

draw.inr <- function(inr.list,high=3.0,low=2.0,label1="",label2="",label3="")
{
    if (length(inr.list$inr)==1) break
    tempnum <- length(inr.list$inr)
    tempdays <- inr.list$day - inr.list$day[1]
    
    plot(tempdays,inr.list$inr,xlab="Day",ylab="INR",ylim=c(0,6),type="l")
    points(tempdays,inr.list$inr,pch=19,cex=1.4)
    abline(h=high,lty=2,lwd=2)
    abline(h=low,lty=2,lwd=2)
    text(3,5,labels=label1,pos=4)
    text(3,4.6,labels=label2,pos=4)
    text(3,5.4,labels=label3,pos=4)
}


# function split.date
# 
# Split dates into the components month, day, and year.

split.date <- function(x,ymd=FALSE,return.year=FALSE,return.month=FALSE,return.day=FALSE,year2d=FALSE)
{
    library(survival)
    splitone <- function(x,ymd=ymd,return.year=return.year,return.month=return.month,return.day=return.day)
    {
        temp <- strsplit(as.character(x),"/")
        month <- rep(0,length(x))
        day <- rep(0,length(x))
        year <- rep(0,length(x))
        res <- rep(NA,length(x))
        if (ymd) {
        if (year2d) year <- 2000 + as.numeric(temp[[1]][1]) else year <- as.numeric(temp[[1]][1])
        month <- as.numeric(temp[[1]][2])
        day <- as.numeric(temp[[1]][3])
        }
        else {
        month <- as.numeric(temp[[1]][1])
        day <- as.numeric(temp[[1]][2])
        if (year2d) year <- 2000 + as.numeric(temp[[1]][3]) else year <- as.numeric(temp[[1]][3])
        }
        if (is.na(day) | is.na(month) | is.na(year)) res <- NA
        else res <- as.numeric(mdy.date(month,day,year))
        if (return.year) return(year) else if (return.month) return (month) else if (return.day) return(day) else return(res)
    }
    return(sapply(x,splitone,ymd=ymd,return.year=return.year,return.month=return.month,return.day=return.day))
}





# function time.in.range
# 
# Calculate TTR using linear interpolation, given the start and end times (t1, t2) and the 
# respective INR values y1,y2.  The upper and lower limits are defined by parameters high 
# and low respectively.

time.in.range <- function(t1,t2,y1,y2,high,low)
{
    t.b <- 0; t.in <- 0; t.a <- 0                                                                    # initialize variables
    auc.b <- 0; auc.a <- 0
    auc2.b <- 0; auc2.a <- 0
    m <- (y2-y1)/(t2-t1-1)                                                                           # slope of line joining two points
    int <- y1-m*(t1+1)                                                                               # intercept of line joining two points
    bh <- (int-high)                                                                                 # deviation from intercept to high limit
    bl <- (low-int)                                                                                  # deviation from intercept to low limit
    if (t2-t1 <= 1) return(c(NA,NA,NA))                                                              # can not proceed

    if (m==0) {                                                                                      # separate case for slope zero
        if (y1 <= high & y1 >= low) t.in <- (t2-t1-1)
        else if (y1 > high) {
            t.a <- (t2-t1-1)
            auc.a <- (y1-high)*t.a
        }
        else if (y1 < low) {
            t.b <- (t2-t1-1)
            auc.b <- (low-y1)*t.b
        }
    }
    else {
    tup <- ((high-int)/m)+t1+1                                                                       # calculate times of intersection with upper and lower bounds
    tdown <- ((low-int)/m)+t1+1

    if (y1 > high & y2 > high) {                                                                     # Case 1: both y1 and y2 above range
        t.a <- (t2-t1-1)                                                                             # entire time above
        auc.a <- (mean(c(y1,y2))-high)*t.a                                                           # trapezoid above
        auc2.a <- (m^2/3)*(t2^3-(t1+1)^3) + m*bh*(t2^2-(t1+1)^2) + bh^2*t.a                          # integrate quadratic over trapezoid
    } 
    else if (y1 < low & y2 < low) {                                                                  # Case 2: both y1 and y2 below range
        t.b <- (t2-t1-1)                                                                             # entire time below
        auc.b <- (low-mean(c(y1,y2)))*t.b                                                            # trapezoid below
        auc2.b <- (m^2/3)*(t2^3-(t1+1)^3) - m*bl*(t2^2-(t1+1)^2) + bl^2*t.b                          # integrate quadratic over trapezoid
    }
    else if ((y1 <= high & y2 <= high & y1 >= low & y2 >= low)) {                                    # Case 3: both y1 and y2 in range
        t.in <- (t2-t1-1)                                                                            # entire time in range
        }
    else if (y1 <= high & y1 >= low & y2 > high) {                                                   # Case 4: y1 in range, y2 above range
        t.in <- tup-t1-1                                                                             # time in range
        t.a <- t2-tup                                                                                # time above
        auc.a <- (t2-tup)*(y2-high)/2                                                                # triangle above
        auc2.a <- (t2-tup)^3*m^2/3                                                                   # integrate quadratic over triangle
        }
    else if (y1 > high & y2 <= high & y2 >= low) {                                                   # Case 5: y1 above range, y2 in range
        t.a <- tup-t1-1                                                                              # time above
        t.in <- t2-tup                                                                               # time in range
        auc.a <- (tup-t1-1)*(y1-high)/2                                                              # triangle above
        auc2.a <- (tup-t1-1)^3*m^2/3                                                                 # integrate quadratic over triangle
        }
    else if (y1 <= high & y1 >= low & y2 < low) {                                                    # Case 6: y1 in range, y2 below range
        t.in <- tdown-t1-1                                                                           # time in range
        t.b <- t2-tdown                                                                              # time below
        auc.b <- (t2-tdown)*(low-y2)/2                                                               # triangle below
        auc2.b <- (t2-tdown)^3*m^2/3                                                                 # integrate quadratic over triangle
        }
    else if (y1 < low & y2 >= low & y2 <= high) {                                                    # Case 7: y1 below range, y2 in range
        t.b <- tdown-t1-1                                                                            # time below
        t.in <- t2-tdown                                                                             # time in range
        auc.b <- (tdown-t1-1)*(low-y1)/2                                                             # triangle below
        auc2.b <- (tdown-t1-1)^3*m^2/3                                                               # integrate quadratic over triangle
        }
    else if (y1 < low & y2 > high) {                                                                 # Case 8: y1 below range, y2 above range
        t.b <- tdown-t1-1                                                                            # time below
        t.in <- tup-tdown                                                                            # time in range
        t.a <- t2-tup                                                                                # time above
        auc.a <- (t2-tup)*(y2-high)/2                                                                # triangle above
        auc.b <- (tdown-t1-1)*(low-y1)/2                                                             # triangle below
        auc2.a <- (t2-tup)^3*m^2/3                                                                   # integrate quadratic over triangle above
        auc2.b <- (tdown-t1-1)^3*m^2/3                                                               # integrate quadratic over triangle below
        }
    else if (y1 > high & y2 < low) {                                                                 # Case 9: y1 above range, y2 below range
        t.a <- tdown-t1-1                                                                            # time above
        t.in <- tdown-tup                                                                            # time in range
        t.b <- t2-tup                                                                                # time below
        auc.a <- (tup-t1-1)*(y1-high)/2                                                              # triangle above
        auc.b <- (t2-tdown)*(low-y2)/2                                                               # triangle below
        auc2.a <- (tup-t1-1)^3*m^2/3                                                                 # integrate quadratic over triangle above
        auc2.b <- (t2-tdown)^3*m^2/3                                                                 # integrate quadratic over triangle below
        }
    }
    return(c(t.b,t.in,t.a,auc.b,auc.a,auc2.b,auc2.a))
}




# function calc.tir
# 
# Calculate TTR from a series of INR measurements.  The parameter 'inr.list' should be a data frame with
# the following format:
#      id     day     inr
#       1      7      2.0
#       1      13     2.5
#       1      18     3.1
#       1      47     2.8

calc.tir <- function(inr.list,lowrange=2.0,highrange=3.0)
{
    inr.list <- inr.list[order(inr.list$id,inr.list$day),]

    idlist <- unique(inr.list$id)
    numsub <- length(idlist)
    resmat <- matrix(NA,numsub,14)
    
    for (i in 1:numsub)
    {
        this.dat <- inr.list[inr.list$id == idlist[i],]
        if (length(this.dat$inr)==1 | sum(is.na(this.dat$day)>0))
        {
            resmat[i,] <- c(idlist[i],0,0,0,0,0,0,NA,NA,NA,NA,NA,NA,NA)
            next
        }
        tempnum <- length(this.dat$inr)
        temp.below <- sum(this.dat$inr < lowrange)
        temp.in <- sum(this.dat$inr <= highrange & this.dat$inr >= lowrange)
        temp.above <- sum(this.dat$inr > highrange)
        temp.gap <- 0
        temp.auca <- 0
        temp.aucb <- 0
        temp.auc2a <- 0
        temp.auc2b <- 0

        for (j in 2:tempnum)
        {
            lastdate <- this.dat$day[j] - this.dat$day[j-1]
            if (lastdate > 56 | lastdate==0 | is.na(this.dat$inr[j-1]) | is.na(this.dat$inr[j])) temp.gap <- temp.gap + (lastdate-1)
            else
            {
                tempres <- time.in.range(0,lastdate,this.dat$inr[j-1],this.dat$inr[j],highrange,lowrange)
                temp.below <- temp.below + round(tempres[1],1)
                temp.in <- temp.in + round(tempres[2],1)
                temp.above <- temp.above + round(tempres[3],1)
                temp.aucb <- temp.aucb + round(tempres[4],1)
                temp.auca <- temp.auca + round(tempres[5],1)
                temp.auc2b <- temp.auc2b + round(tempres[6],1)
                temp.auc2a <- temp.auc2a + round(tempres[7],1)
            }
        }
        total.time <- temp.below+temp.in+temp.above+temp.gap
        time.nogaps <- total.time-temp.gap
        resmat[i,] <- c(idlist[i],total.time,temp.gap,time.nogaps,temp.below,temp.in,temp.above,round(temp.in/time.nogaps,3),
                         round(temp.below/time.nogaps,3),round(temp.above/time.nogaps,3),auc.below=temp.aucb,temp.auca,
                         temp.auc2b,temp.auc2a)
    }
    
    resmat <- data.frame(resmat)
    names(resmat) <- c("id","total.time","gaps","time.nogaps","time.below","time.in","time.above","tir","tir.below","tir.above","auc.below","auc.above","auc2.below","auc2.above")
    return(resmat)

#    return(list(total.time=total.time,
#                gaps=temp.gap,
#                time.nogaps=time.nogaps,
#                time.below=temp.below,
#                time.in=temp.in,
#                time.above=temp.above,
#                tir=round(temp.in/time.nogaps,3),
#                tir.below=round(temp.below/time.nogaps,3),
#                tir.above=round(temp.above/time.nogaps,3),
#                auc.below=temp.aucb,
#                auc.above=temp.auca,
#                auc2.below=temp.auc2b,
#                auc2.above=temp.auc2a))

}



# function calc.sigma
#
# Calculate Fihn's sigma from a series of INR measurements.  The parameter this.inr should be
# a data frame of INR values in identical format to that described above.  Past data suggest this
# variable typically follows a log-normal distribution.

calc.sigma <- function(inr.list)
{
    num <- length(inr.list$inr)
    if (num < 2) return(NA)
    temp <- 0
    for (j in 2:num)
    {
        lastdate <- inr.list$day[j] - inr.list$day[j-1]
        if (lastdate < 56) temp <- temp + ((inr.list$inr[j] - inr.list$inr[j-1])^2 / (lastdate/7))
    }
    sigma <- sqrt(temp/(num - 1))
    return(replace(sigma,sigma==0 | sigma==Inf,NA))
}
