library(latChanNet)
countList = data.frame(i = c(1,1,2,3, 1), 
                       j = c(2, 3, 3, 4, 5), 
                       cnt = c(1,1, 1, 1, 1))


bkn_mod = makeBKN(countList, 2)
bkn_mod$get_theta()
bkn_mod$llk()
emBKN(bkn_mod, 100)
bkn_mod$llk()
bkn_mod = makeBKN(countList, 2)
emBKN(bkn_mod, 1, type = 2)
bkn_mod$llk()
bkn_mod$get_theta()
