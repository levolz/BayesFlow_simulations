# get priors out of `corey` object in "corey_noError_hier.RData"
priors = list()
for(i in seq_along(corey)){
    for (j in names(corey[[i]]$p.prior)){
        if (i == 1) priors[[j]] = unlist(corey[[i]]$p.prior[[j]])
        else priors[[j]] = rbind(priors[[j]], unlist(corey[[i]]$p.prior[[j]]))
    }
}

# prior means from `corey` fits