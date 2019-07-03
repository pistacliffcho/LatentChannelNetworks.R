library(latChanNet)
countList = data.frame(i = c(1, 1, 2, 3, 1), 
                       j = c(2, 3, 3, 4, 5), 
                       cnt = c(1,1, 1, 1, 1))


bkn_mod = makeBKN(countList, 2)
bkn_mod$get_theta()
bkn_mod$llk()
emBKN(bkn_mod, 1000)
bkn_mod$llk()

missing_list = data.frame(i = c(1, 2), 
                          j = c(4, 4))
bkn_mod2 = makeBKN(countList, 2, missing_list)
emBKN(bkn_mod2, 1000)
bkn_mod2$llk()
