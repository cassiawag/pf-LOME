library("R6")
source("PfLOME_Pathogen.R")
source("PfLOME_Human.R")
source("Rx.R")


dt = 1
nptypes = c(3,4,6)
pfid = 0
someGuy = Human$new(1)

bites = 0
tt = 0
while(tt<365*2){
  bite = rgeom(1,100/365)
  # rgeom is also a memoryless distribution, but discrete unlike continuous exponentioal
  # rate used her is 100 bites per year on average
  bites = c(bites,bite)
  tt = cumsum(bites)[length(bites)]
}
bites = unique(cumsum(bites))
#bites = unique(sort(make.bites(70, 10, 1, 5, wt=wt, trend = .05)))
moi = 1+rnbinom(length(bites), mu=3, size = .3)

#dt = 7
tFinal = 2*365 #10
t = 1

while(t < tFinal){
  someGuy$updateHuman(t,dt)
  s = t
  while((s >= t) & s < t+dt){
    if(s %in% bites){
      k = which(bites==s)
      for(i in 1:moi[k]){
        pfid = pfid + 1
        pf = Pf$new(pfid,nptypes) 
        someGuy$infectHuman(s,pf)
      }
    }
    s = s+1
  }
  if(someGuy$get_Fever()>0){
    p = rbinom(1,1,.03)
    if(p == 1){
      someGuy$Treat(t,1)
    }
  }
  t = t+dt
}
######################### plotting functions #############################

tFinal = ifelse(dt%%2==1,tFinal-1,tFinal)
t = seq(1,tFinal,dt)/365
plot(t,someGuy$get_history()$Ptot,type="l",
     ylim=c(-5,11),xlim=c(0,2),xlab='years',ylab='log10 iRBC')
lines(t,someGuy$get_history()$Gtot,lty=2)
lines(t,someGuy$get_history()$Fever)
lines(t,someGuy$get_history()$GenImm,type="l") # get error about x & y lengths differ. Tried to fix by removing -3 & 2*
lines(t,someGuy$get_history()$PD,col='purple')
abline(h=c(-1,-3),lty=2)
lines(t,someGuy$get_history()$PfMOI/max(someGuy$get_history()$PfMOI)*2-5)
abline(h=-5,lty=2)
text(.1,-3.5, paste("MOI, max=", max(someGuy$get_history()$PfMOI)), col = "blue", pos = 4)

plot(t,someGuy$get_history()$GenImm,type="l",xlab='years',ylab='% of max strength of immunity')
for(i in 1:someGuy$get_immuneState()$get_nBSImmCounters()){
  lines(c(0,t),someGuy$get_history()$BSImm[[i]],lty=2)
}

plot(1:length(someGuy$get_history()$PfMOI),someGuy$get_history()$PfMOI,type="l",
     xlab='days',ylab='MOI')

##RBC
plot(1:length(someGuy$get_history()$RBC),someGuy$get_history()$RBC,type="l",xlab ='days')
##HRP2
plot(1:length(someGuy$get_history()$HRP2),someGuy$get_history()$HRP2,type="l",xlab='days')
##pLDH
plot(1:length(someGuy$get_history()$pLDH),someGuy$get_history()$pLDH,type="l",xlab='days')

