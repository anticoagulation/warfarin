
##########
#
# Functions to calculate various measures of INR control, including
# Time in Therapeutic Range (TTR) following the linear interpolation 
# method of Rosendaal.
#
# Project started 2014-03-05 by Al Ozonoff, Boston Children's Hospital
# 2014-03-17 made changes to accept multiple patient id's
# 2014-05-02 further changes to accept missing or consecutive INR values
# 2014-08-25 fixed problems calculating intercept, Cases 4/5/6/7 of calc.tir
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

time.in.range <- function(t1,t2,y1,y2,high=3.0,low=2.0)
{
    t.b <- 0; t.in <- 0; t.a <- 0                                                                    # initialize variables
    auc.b <- 0; auc.a <- 0
    auc2.b <- 0; auc2.a <- 0
    m <- (y2-y1)/(t2-t1-1)                                                                           # slope of line joining two points
    int <- y1-m*t1                                                                                   # intercept of line joining two points
    bh <- (int-high)                                                                                 # deviation from intercept to high limit
    bl <- (low-int)                                                                                  # deviation from intercept to low limit
    if (t2-t1 < 1) return(rep(NA,7))                                                                 # two INR values less than one day apart - something is wrong
    if (t2-t1 == 1) return(rep(0,7))                                                                 # two INR values one day apart - recalculate

    if (m==0) {                                                                                      # separate case for slope zero
        if (y1 <= high & y1 >= low) t.in <- (t2-t1-1)
        else if (y1 > high) {
            t.a <- (t2-t1-1)
            auc.a <- (y1-high)*t.a
            auc2.a <- (y1-high)^2*t.a
        }
        else if (y1 < low) {
            t.b <- (t2-t1-1)
            auc.b <- (low-y1)*t.b
            auc2.b <- (low-y1)^2*t.b
        }
    }
    else {
    tup <- ((high-int)/m)                                                                        # calculate times of intersection with upper and lower bounds
    tdown <- ((low-int)/m)

    if (y1 > high & y2 > high) {                                                                     # Case 1: both y1 and y2 above range
        t.a <- (t2-t1-1)                                                                             # entire time above
        auc.a <- (mean(c(y1,y2))-high)*t.a                                                           # trapezoid above
        auc2.a <- ((m^2/3)*t.a^3) + (y1-high)^2*t.a                                                  # integrate quadratic over trapezoid
    }
    else if (y1 < low & y2 < low) {                                                                  # Case 2: both y1 and y2 below range
        t.b <- (t2-t1-1)                                                                             # entire time below
        auc.b <- (low-mean(c(y1,y2)))*t.b                                                            # trapezoid below
        auc2.b <- ((m^2/3)*t.b^3) + (low-y1)^2*t.b                                                   # integrate quadratic over trapezoid
    }
    else if ((y1 <= high & y2 <= high & y1 >= low & y2 >= low)) {                                    # Case 3: both y1 and y2 in range
        t.in <- (t2-t1-1)                                                                            # entire time in range
        }
    else if (y1 <= high & y1 >= low & y2 > high) {                                                   # Case 4: y1 in range, y2 above range
        t.in <- tup-t1                                                                               # time in range
        t.a <- t2-tup-1                                                                              # time above
        auc.a <- (t2-tup-1)*(y2-high)/2                                                              # triangle above
        auc2.a <- (t2-tup-1)^3*m^2/3                                                                 # integrate quadratic over triangle
        }
    else if (y1 > high & y2 <= high & y2 >= low) {                                                   # Case 5: y1 above range, y2 in range
        t.in <- t2-tup-1                                                                             # time in range
        t.a <- tup-t1                                                                                # time above
        auc.a <- (tup-t1)*(y1-high)/2                                                              # triangle above
        auc2.a <- (tup-t1)^3*m^2/3                                                                 # integrate quadratic over triangle
        }
    else if (y1 <= high & y1 >= low & y2 < low) {                                                    # Case 6: y1 in range, y2 below range
        t.b <- t2-tdown-1                                                                            # time below
        auc.b <- (t2-tdown-1)*(low-y2)/2                                                             # triangle below
        auc2.b <- (t2-tdown-1)^3*m^2/3                                                               # integrate quadratic over triangle
        }
    else if (y1 < low & y2 >= low & y2 <= high) {                                                    # Case 7: y1 below range, y2 in range
        t.b <- tdown-t1                                                                              # time below
        auc.b <- (tdown-t1)*(low-y1)/2                                                               # triangle below
        auc2.b <- (tdown-t1)^3*m^2/3                                                                 # integrate quadratic over triangle
        }
    else if (y1 < low & y2 > high) {                                                                 # Case 8: y1 below range, y2 above range
        t.b <- tdown-t1                                                                              # time below
        t.in <- tup-tdown                                                                            # time in range
        t.a <- t2-tup-1                                                                              # time above
        auc.a <- (t2-tup-1)*(y2-high)/2                                                              # triangle above
        auc.b <- (tdown-t1)*(low-y1)/2                                                               # triangle below
        auc2.a <- (t2-tup-1)^3*m^2/3                                                                 # integrate quadratic over triangle above
        auc2.b <- (tdown-t1)^3*m^2/3                                                                 # integrate quadratic over triangle below
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
    inr.list <- inr.list[!is.na(inr.list$inr),]
    inr.list <- inr.list[order(inr.list$id,inr.list$day),]

    idlist <- unique(inr.list$id)
    numsub <- length(idlist)
    resmat <- matrix(NA,numsub,14)
    
	total.time.in.allpts = 0 		#Bob: initialized variables to collect cumulative data across patients 03-19-2014
	total.time.nogaps.allpts = 0
	total.time.allpts = 0

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
        temp.auca <- sum(as.numeric(this.dat$inr > highrange)*sapply((this.dat$inr-highrange),max,0))
        temp.aucb <- sum(as.numeric(this.dat$inr < lowrange)*sapply((lowrange-this.dat$inr),max,0))
        temp.auc2a <- sum(as.numeric(this.dat$inr > highrange)*(sapply((this.dat$inr-highrange),max,0))^2)
        temp.auc2b <- sum(as.numeric(this.dat$inr < lowrange)*(sapply((lowrange-this.dat$inr),max,0))^2)

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
		total.time.in.allpts = total.time.in.allpts + temp.in				#Bob added to increment variables to collect cumulative data across patients 03-19-2014
		total.time.nogaps.allpts = total.time.nogaps.allpts + time.nogaps
		total.time.allpts = total.time.allpts + total.time
    }
    
	#Return the totals of the variables that collected cumulative data across patients 03-19-2014
	return(list(tir=round(total.time.in.allpts / total.time.nogaps.allpts,3),total.time.nogaps.allpts = total.time.nogaps.allpts, total.time.allpts = total.time.allpts, gaps.fraction = round(1 - (total.time.nogaps.allpts / total.time.allpts),3), number.subjects = numsub))
    #resmat <- data.frame(resmat)
    #names(resmat) <- c("id","total.time","gaps","time.nogaps","time.below","time.in","time.above","tir","tir.below","tir.above","auc.below","auc.above","auc2.below","auc2.above")
    #allttr <- round(sum(resmat$time.in,na.rm=T)/sum(resmat$time.nogaps,na.rm=T),3)
    #numsub <- sum(!is.na(resmat$tir))
    #alltime <- sum(resmat$time.nogaps,na.rm=T)
    #return(resmat)

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

# function launch
# 
# Accept input from web form and call calc.tir 
launch <- function(content, lowrange, highrange) {
temp <- content
# http://stat.ethz.ch/R-manual/R-devel/library/base/html/regex.html
#temp <- gsub('\n', '', fixed = TRUE, temp, perl = TRUE)
#temp <- gsub("\\s+$", "", temp, perl = TRUE) #Removing trailing whitespace
temp <- gsub(",+$", "", temp, perl = TRUE) #Remove trailing comma if accidentally added by user online at webform
temp <- gsub("\t", ' ', fixed = TRUE, temp)
temp <- gsub(',', '","', fixed = TRUE, temp)
temp <- paste('"',temp,'"',sep = '')
temp <- paste('Mymatrix <- matrix(c(',temp,'), ncol=3, byrow=TRUE, dimnames = list(NULL, c("id","day", "inr")))',sep = '')
x <- eval(parse(file = "", n = NULL, text = temp))
inr.list <- data.frame (x)
inr.list$id  <- as.numeric(as.character(str_trim(inr.list$id)))
inr.list$day <- as.numeric(as.character(str_trim(inr.list$day)))
inr.list$inr <- as.numeric(as.character(str_trim(inr.list$inr)))
lowrange <- as.numeric(as.character(str_trim(lowrange)))
highrange <- as.numeric(as.character(str_trim(highrange)))
#Shebang
result <- calc.tir (inr.list, lowrange, highrange)
svgtext = "<?xml version=\"1.0\" encoding=\"UTF-8\" standalone=\"no\"?><!-- Created with Inkscape (http://www.inkscape.org/) --><svg   xmlns:dc=\"http://purl.org/dc/elements/1.1/\"   xmlns:cc=\"http://creativecommons.org/ns#\"   xmlns:rdf=\"http://www.w3.org/1999/02/22-rdf-syntax-ns#\"   xmlns:svg=\"http://www.w3.org/2000/svg\"   xmlns=\"http://www.w3.org/2000/svg\"   version=\"1.2\"   width=\"650\"   height=\"600\"   id=\"svg\">  <defs     id=\"defs11\" />  <metadata     id=\"metadata2991\">    <rdf:RDF>      <cc:Work         rdf:about=\"\">        <dc:format>image/svg+xml</dc:format>        <dc:type           rdf:resource=\"http://purl.org/dc/dcmitype/StillImage\" />        <dc:title></dc:title>      </cc:Work>    </rdf:RDF>  </metadata>  <rect     width=\"200\"     height=\"200\"     x=\"200\"     y=\"0\"     id=\"high.color\"     style=\"fill:#00ff00\" />  <rect     width=\"200\"     height=\"200\"     x=\"200\"     y=\"200\"     id=\"medium.color\"     style=\"fill:#ffcc00\" />  <rect     width=\"200\"     height=\"200\"     x=\"200\"     y=\"400\"     id=\"low.color\"     style=\"fill:#ff8080\" />  <text     x=\"5\"     y=\"202\"     id=\"text.best\"     xml:space=\"preserve\"><tspan       x=\"5\"       y=\"202\"       id=\"tspan3112\"       style=\"font-size:16px;font-weight:bold\">Best practice</tspan><tspan       x=\"5\"       y=\"222\"       id=\"tspan3125\"       style=\"font-size:16px\">(randomized trials): 75%</tspan></text>  <polygon     points=\"233.17773,364.71875 0,182.35938 233.17773,0 233.17773,141.73242 902.25049,141.73242 902.25049,222.98633 902.25049,222.98633 233.17773,222.98633 \"     transform=\"matrix(-0.07575744,0,0,-0.05197441,201.30134,208.71875)\"     id=\"path-arrow-best\" />  <text     x=\"5\"     y=\"412\"     id=\"text.typical\"     xml:space=\"preserve\"><tspan       x=\"5\"       y=\"437\"       id=\"tspan3112-4\"       style=\"font-size:16px;font-weight:bold\">Typical practice</tspan><tspan       x=\"5\"       y=\"457\"       id=\"tspan3125-5\"       style=\"font-size:16px\">(community setting): 58%</tspan></text>  <polygon     points=\"902.25049,141.73242 902.25049,222.98633 902.25049,222.98633 233.17773,222.98633 233.17773,364.71875 0,182.35938 233.17773,0 233.17773,141.73242 \"     transform=\"matrix(-0.07575744,0,0,-0.05197441,200.34559,439.78311)\"     id=\"path-arrow-typical\" />  <text     x=\"250\"     y=\"19.476227\"     id=\"TTR_high\"     style=\"font-size:18px;font-weight:bold\">TTR &gt; 75%</text>  <text     x=\"210\"     y=\"76.816208\"     id=\"text-high-1\"     xml:space=\"preserve\"     style=\"font-size:14px\">Per 100 patients-years:</text>  <text     x=\"210\"     y=\"94.823372\"     id=\"text-high-2\"     xml:space=\"preserve\"     style=\"font-size:14px\">Stroke 1%</text>  <text     x=\"210\"     y=\"113.9599\"     id=\"text-high-3\"     xml:space=\"preserve\"     style=\"font-size:14px\">Major bleed 1%</text>  <text     x=\"210\"     y=\"134.53052\"     id=\"text-high-4\"     xml:space=\"preserve\"     style=\"font-size:14px\">Death 2%</text>  <text     x=\"230\"     y=\"221.47623\"     id=\"TTR_medium\"     style=\"font-size:18px;font-weight:bold\">TTR 60% to 75%</text>  <text     x=\"210\"     y=\"270.91919\"     id=\"text-medium-1\"     xml:space=\"preserve\"     style=\"font-size:14px\">Per 100 patients-years:</text>  <text     x=\"210\"     y=\"289.8252\"     id=\"text-medium-2\"     xml:space=\"preserve\"     style=\"font-size:14px\">Stroke 1%</text>  <text     x=\"210\"     y=\"308.96173\"     id=\"text-medium-3\"     xml:space=\"preserve\"     style=\"font-size:14px\">Major bleed 2%</text>  <text     x=\"210\"     y=\"328.63348\"     id=\"text-medium-4\"     xml:space=\"preserve\"     style=\"font-size:14px\">Death 2%</text>  <text     x=\"250\"     y=\"420.24124\"     id=\"TTR_low\"     style=\"font-size:18px;font-weight:bold\">TTR &lt; 60%</text>  <text     x=\"210\"     y=\"476.11118\"     id=\"text-low-1\"     xml:space=\"preserve\"     style=\"font-size:14px\">Per 100 patients-years:</text>  <text     x=\"210\"     y=\"492.80011\"     id=\"text-low-2\"     xml:space=\"preserve\"     style=\"font-size:14px\">Stroke 2%</text>  <text     x=\"210\"     y=\"510.70737\"     id=\"text-low-3\"     xml:space=\"preserve\"     style=\"font-size:14px\">Major bleed 4%</text>  <text     x=\"210\"     y=\"526.38239\"     id=\"text-low-4\"     xml:space=\"preserve\"     style=\"font-size:14px\">Death 4%</text>  <polygon     points=\"0,182.35938 233.17773,0 233.17773,141.73242 902.25049,141.73242 902.25049,222.98633 902.25049,222.98633 233.17773,222.98633 233.17773,364.71875 \"     transform=\"matrix(0.06375505,0.04092045,0.02807403,-0.04374001,396.67506,433.02394)\"     id=\"path-arrow-threshold\" />  <text     x=\"430\"     y=\"480\"     id=\"flowPara4477\"     style=\"font-size:14px;font-weight:bold\">TTR = 58%</text>  <text     x=\"430\"     y=\"500\"     id=\"flowPara4479\">Threshold modeled for</text>  <text     x=\"430\"     y=\"515\"     id=\"flowPara3019\">benefit over ASA/clopidogrel in</text>  <text     x=\"430\"     y=\"530\"     id=\"flowPara3025\">atrial fibrillation</text>  <rect     width=\"225\"     height=\"80\"     x=\"415\"     y=\"460\"     id=\"rect_threshold\"     style=\"fill:#ffffff;fill-opacity:0;stroke:#000000;stroke-width:2.25650907;stroke-miterlimit:4;stroke-opacity:1;stroke-dasharray:none\" />  <text     x=\"405\"     y=\"207.1962\"     id=\"text4755\"     xml:space=\"preserve\"     style=\"font-size:24px;font-style:normal;font-weight:normal;line-height:125%;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;font-family:Sans\">75%</text>  <text     x=\"405\"     y=\"408.36755\"     id=\"text4755-4\"     xml:space=\"preserve\"     style=\"font-size:24px;font-style:normal;font-weight:normal;line-height:125%;letter-spacing:0px;word-spacing:0px;fill:#000000;fill-opacity:1;stroke:none;font-family:Sans\">60%</text>  <text     x=\"50\"     y=\"100\"     id=\"result-high\"     xml:space=\"preserve\"     style=\"font-size:40px;font-style:normal;font-variant:normal;font-weight:bold;font-stretch:normal;text-align:start;line-height:125%;writing-mode:lr-tb;text-anchor:start;fill:#00ff00;font-family:Sans;-inkscape-font-specification:Sans Bold\">HH%</text>  <text     x=\"50\"     y=\"310\"     id=\"result-medium\"     xml:space=\"preserve\"     style=\"font-size:40px;font-style:normal;font-variant:normal;font-weight:bold;font-stretch:normal;text-align:start;line-height:125%;writing-mode:lr-tb;text-anchor:start;fill:#ffcc00;font-family:Sans;-inkscape-font-specification:Sans Bold\">MM%</text>  <text     x=\"50\"     y=\"510\"     id=\"result-low\"     xml:space=\"preserve\"     style=\"font-size:40px;font-style:normal;font-variant:normal;font-weight:bold;font-stretch:normal;text-align:start;line-height:125%;writing-mode:lr-tb;text-anchor:start;fill:#ff8080;font-family:Sans;-inkscape-font-specification:Sans Bold\">LL%</text>  <polygon     points=\"0,182.35938 233.17773,0 233.17773,141.73242 902.25049,141.73242 902.25049,222.98633 902.25049,222.98633 233.17773,222.98633 233.17773,364.71875 \"     transform=\"matrix(0.07346576,-0.02419435,0.02226995,0.04433791,397.51679,327.41265)\"     id=\"path-arrow-threshold-1\" />  <text     x=\"430\"     y=\"250\"     id=\"flowPara4477-0\"     style=\"font-size:14px;font-weight:bold\">TTR = 65%</text>  <text     x=\"430\"     y=\"270\"     id=\"flowPara4479-7\">Threshold modeled for</text>  <text     x=\"430\"     y=\"285\"     id=\"flowPara4479-7\">benefit over DOACs in</text>  <text     x=\"430\"     y=\"300\"     id=\"flowPara4479-7\">atrial fibrillation</text>  <rect     width=\"225\"     height=\"80\"     x=\"415\"     y=\"230\"     id=\"rect_threshold-4\"     style=\"fill:#ffffff;fill-opacity:0;stroke:#000000;stroke-width:2.25650907;stroke-miterlimit:4;stroke-opacity:1;stroke-dasharray:none\" /></svg>"
if(100*result$tir >75){
svgtext <- sub("HH", sprintf("%.1f",100*result$tir),svgtext)
svgtext <- sub("MM%","",svgtext)
svgtext <- sub("LL%","",svgtext)
} else if(100*result$tir <60){
svgtext <- sub("HH%","",svgtext)
svgtext <- sub("MM%","",svgtext)
svgtext <- sub("LL",sprintf("%.1f",100*result$tir),svgtext)
} else {
svgtext <- sub("HH%","",svgtext)
svgtext <- sub("MM", sprintf("%.1f",100*result$tir),svgtext)
svgtext <- sub("LL%", "",svgtext)
}
msg <- paste("<div>Time in therapeutic range: ", sprintf("%.1f",100*result$tir), "%</div>", sep="")
msg <- paste(svgtext,msg,"<div>Total number of days without gaps: ", sprintf("%.1f",result$total.time.nogaps.allpts), ". Total time: ", sprintf("%.1f",result$total.time.allpts), ". Gaps fraction: ", sprintf("%.1f",100*result$gaps.fraction), "%</div>", sep="")
msg <- paste(msg,"<div>Total number of subjects: ", result$number.subjects, "</div>", sep="")

list(message = msg)
}
