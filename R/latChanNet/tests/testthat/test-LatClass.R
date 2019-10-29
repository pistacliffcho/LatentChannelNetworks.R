data("email_data")

set.seed(1)

# Testing basic model w/out metadata
mod = makeLatentModel(email_data$edgeList, 10)
llk_before = mod$llk()
alg_res <- mod$fit(fast_em = T)
llk_after = mod$llk()
llk_greater = llk_after > llk_before
test_that("Base LCN EM", 
          {expect_true(llk_greater)})

# Testing predictions
i = 1:2; j = 3:4
cross_preds = predict(mod, i, j, "cross")
row1_preds = predict(mod, i[1], j, "pairs")
row2_preds = predict(mod, i[2], j)
comb_preds = rbind(row1_preds, row2_preds)
all_eq = all(cross_preds == comb_preds)
test_that("Cross prediction == row by row predictions", {
  expect_true(all_eq)})


# Testing metadata
mod = makeLatentModel(email_data$edgeList, 10, metadata = email_data$meta)
alg_res <- mod$fit(fast_em = T)
