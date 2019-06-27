library(latChanNet)
countList = data.frame(i = c(1,1,2,3), 
                       j = c(2, 3, 3, 4), 
                       cnt = c(1,1, 2, 3))


bkn_mod = makeBKN(countList, 2)
bkn_mod$get_theta()
bkn_mod$llk()
for(i in 1:1000)
  bkn_mod$one_em()
bkn_mod$llk()
