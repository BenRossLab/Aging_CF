

improve_plot <- function() {
  #####
  # Add to any ggplot2 plot to make it look better
  #####
  
  library('ggplot2')
  t <- theme(axis.text = element_text(color = "black"),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.x = element_blank(),
        panel.grid.minor.y = element_blank(),
        panel.grid.major.y = element_line(linetype="dashed"))
  return(t)
}



plot_species_over_time <- function(df, sp) {
    p <- ggplot(df, aes(y = get(sp), x = age)) + geom_point(pch=21, size=3,aes(fill = get(sp) > 0)) + theme_bw(16) + labs(fill = "Detectable") + 
                                       xlab("Age (months)") + theme(axis.text = element_text(color="black")) + ylab(sprintf("%s (abundance)", str_replace(sp, "s__","")))
    
    cor.val <- cor(pull(df_model_sub, get(sp)), pull(df_model_sub,age), method = "spearman")
    return(p + ggtitle(sprintf("spearman correlation %.2f", cor.val)))
 
}

plot_species_discrete_break <- function(df, sp, levels_names, bin_width=0.25, type = "fill", fill_color = "steelblue") {

    t <- improve_plot()

    df$old <- factor(df$old)
    levels(df$old) <- levels_names

    if (type == "facet") {
        # https://stackoverflow.com/questions/52690318/percentage-histogram-with-facet-wrap
        # Justification for how to get percentage of data
        p <- ggplot(df, aes(get(sp), group=old)) + geom_histogram(aes(y = stat(density) * bin_width), binwidth=bin_width, fill=fill_color, color = "black") + facet_grid(old ~ .) + theme_bw(16) + xlab(paste0(sp, " (abundance)")) + ylab("Percentage of data") + scale_y_continuous(labels = scales::percent) + t
    } else if (type == "density") {
        p <- ggplot(df, aes(get(sp), group=old)) + geom_density(aes(fill=old), alpha=0.5) + theme_bw(16) + xlab(paste0(sp, " (abundance)")) + ylab("Density of data") + t

    } else {
        p <- ggplot(df, aes(get(sp), group=old)) + geom_histogram(aes(y = c(..count..[..group..==1]/sum(..count..[..group..==1]),
                         ..count..[..group..==2]/sum(..count..[..group..==2])), fill=old), binwidth=bin_width, position=position_dodge(), color = "black") + theme_bw(16) + xlab(paste0(sp, " (abundance)")) + ylab("Percentage of data") + scale_y_continuous(labels = scales::percent) + scale_fill_manual(values = c("steelblue","#c50606"), name="age range (months)") +  t
    }
    return(p)
}

get_importance_df <- function(model, type=2) {
    df_importance <- data.frame(RRF::importance(model, type=type))
    df_importance$species <- rownames(df_importance)
    df_importance <- as_tibble(df_importance)
    if (type == 2) {
        df_importance <- arrange(df_importance, desc(IncNodePurity))
    } else {
        df_importance <- arrange(df_importance, desc(X.IncMSE))
    }
    return(df_importance)
}

predict_in_cv <- function(df, response_var = "age", n.flds=10,  verbose=FALSE, ...) {
    library('splitTools')
    library('RRF')

    predictions <- vector(mode = "numeric", length=nrow(df))

    flds <- create_folds(pull(df, response_var), k=n.flds, type = "stratified")
    for (i in 1:n.flds) {
        if (verbose)
            print(sprintf("Currently on fold %i",i))
        idx.train <- flds[[i]]
        idx.test <- 1:nrow(df)
        idx.test <- idx.test[!idx.test%in%idx.train]

        model <- RRF(reformulate(termlabels = '.', response = response_var), df[idx.train,], ...)
        predictions[idx.test] <- predict(model, df[idx.test,])
    }
    return(predictions)
}


build_submodels <- function(df_importance, df_model_sub, n.plot=20, y.var="old", metric ="roc", ntree=100) {
    if (metric == "roc")
        metric.best <- roc.area(df_model_sub[,y.var], predict_in_cv(df_model_sub, y.var))$A
    else if (metric == "cor")
        metric.best <- cor(df_model_sub[,y.var], predict_in_cv(df_model_sub, y.var), method = "spearman")

    df_submodel <- data.frame()
    species.list <- as.character(df_importance$species)
    for (n in 1:n.plot) {
        species.oi <- species.list[1:n]

        df_model_sub_species <- df_model_sub[,c(species.oi, y.var)]

        # model <- RRF(as.formula(sprintf("%s ~ .", y.var)), df_model_sub_species, ntree=ntree)
        # if (metric == "roc")
        #     df_submodel <- bind_rows(df_submodel, data.frame(n=n, roc=roc.area(df_model_sub_species[,y.var], predict(model))$A))
        # else if (metric == "cor")
        #     df_submodel <- bind_rows(df_submodel, data.frame(n=n, cor=cor(df_model_sub_species[,y.var], predict(model), method = "spearman")))
        if (metric == "roc")
            df_submodel <- bind_rows(df_submodel, data.frame(n=n, roc=roc.area(df_model_sub_species[,y.var], predict_in_cv(df_model_sub_species,y.var))$A))
        else if (metric == "cor")
            df_submodel <- bind_rows(df_submodel, data.frame(n=n, cor=cor(df_model_sub_species[,y.var], predict_in_cv(df_model_sub_species,y.var), method = "spearman")))

    }
    df_submodel$species <- species.list[1:n.plot]
    df_submodel$species <- factor(df_submodel$species, levels=df_submodel$species, ordered=TRUE)

    metric.use.ratio <- paste0(metric,".ratio")
    df_submodel[,metric.use.ratio] <- pull(df_submodel,get(metric)) / metric.best

    return(df_submodel)
}


add_roc_title <- function(df_model_sub, sp) {
    vals <- pull(df_model_sub, get(sp))

    roc.val <- roc.area(df_model_sub$old, vals)$A

    p.val <- wilcox.test(x = vals[df_model_sub$old], y = vals[!df_model_sub$old])$p.value
    if (roc.val < 0.5)
        roc.val <- 1.0-roc.val
    return(sprintf("ROC: %.2f, pval: %.0e", roc.val, p.val))
}

create_age_split_df <- function(df_model, low, m, high) {

    df_model_sub <- dplyr::filter(df_model, (age >= low) & (age <= high))
    df_model_sub$old <- df_model_sub$age >= m
    df_model_sub <- dplyr::select(df_model_sub,-age)
    return(df_model_sub)
}



find_critical_species <- function(df_model, method, cutoffs.use = list(c(0,6,12),c(7,13,19)), cutoff.fraction.metric=0.90) {
    source("Lib/Figures.R")

    library('RRF')
    library('tidyverse')
    library('verification')

    t <- improve_plot()

    df_species_use <- NULL


    if (method == "quantitative") {
        response_var <- "age"
        metric <- "cor"
    }
    else {
        response_var <- "old"
        metric <- "roc"
    }
    metric.ratio <- paste0(metric,".ratio")

    for (i in 1:length(cutoffs.use)) {
        cutoffs <- cutoffs.use[[i]]
        low <- cutoffs[1]
        m <- cutoffs[2]
        high <- cutoffs[3]

        df_model_sub <- filter_age_range(df_model, method, low, m, high)
        predictions <- predict_in_cv(df_model_sub, response_var=response_var)
        
        val <- eval_value(predictions, df_model_sub, m)
        print(sprintf("Cutff %.2f: Best CV value is %.2f", m, cor.val))

        model <- RRF(as.formula(sprintf("%s ~ .",response_var)), df_model_sub)
        df_importance <- get_importance_df(model)
        
        df_plot <- build_submodels(df_importance, df_model_sub, y.var = response_var, metric=metric)
        df_plot[,metric.ratio] <- df_plot[,metric] / val
        p <- ggplot(df_plot, aes_string(x = "species", y = metric)) + geom_bar(stat="identity", fill="steelblue", color = "black") + theme_bw(16) + theme(axis.text.x = element_text(angle = 90, vjust = 0.5, hjust=1)) + t
        print(p)

        n.use <- filter(df_plot, get(metric.ratio) >= cutoff.fraction.metric)$n[1]
        species.use <- as.character(df_plot[1:n.use,]$species)
        df_species_use <- bind_rows(df_species_use, data.frame(m=m, species=species.use, metric=df_plot[,metric][1:n.use]))
    }
    return(df_species_use)
}



eval_value <- function(predictions, df_model_sub, method, cor.method = "spearman") {
    require('verification')

    if (method == "quantitative") {
        cor.val <- cor(predictions, df_model_sub$age, method = cor.method)
        return(cor.val)
    } else {
        roc.val <- roc.area(df_model_sub$old, predictions)$A
        return(roc.val)
    }

}


filter_age_range <- function(df_model, method, low, m, high) {
    if (method == "quantitative")
        df_model_sub <- dplyr::filter(df_model, (age >= low) & (age <= high))
    else
        df_model_sub <- create_age_split_df(df_model, low, m, high)
    df_model_sub <- na.omit(df_model_sub)
    return(df_model_sub)
}


cross_dataset_model <- function(df_model1, df_model2, name1, name2, method = "quantitative", cor.method="spearman",ntree=10000) {

    if (method == "quantitative")
        response_var <- "age"
    else
        response_var <- "old"

    print(response_var)
    # Within dataset predictions
    predictions.1 <- predict_in_cv(df_model1, response_var=response_var, ntree=ntree)
    val.1.within<- eval_value(predictions.1, df_model1, method, cor.method)
    
    predictions.2 <- predict_in_cv(df_model2, response_var=response_var, ntree=ntree)
    val.2.within <- eval_value(predictions.2, df_model2, method, cor.method)
    
    # Cross dataset predictions
    model.1 <- RRF(as.formula(sprintf("%s ~ .",response_var)), df_model1)
    model.2 <- RRF(as.formula(sprintf("%s ~ .",response_var)), df_model2)

    predictions.1.cross <- predict(model.2, df_model1)
    predictions.2.cross <- predict(model.1, df_model2)
    val.1.cross <- eval_value(predictions.1.cross, df_model1, method, cor.method)
    val.2.cross <- eval_value(predictions.2.cross, df_model2, method, cor.method)

    df_results <- data.frame(dataset = c(name1, name2, name1, name2), 
                      correlations=c(val.1.within, val.2.within, val.1.cross, val.2.cross),
                      ml = c("self_prediction","self_prediction","cross_prediction","cross_prediction"))

    
    results <- list()
    results$df_results <- df_results
    results$df_plot.within.1 <- data.frame(actual=pull(df_model1, get(response_var)), prediction=predictions.1)
    results$df_plot.cross.1 <- data.frame(actual=pull(df_model1, get(response_var)), prediction=predictions.1.cross)
    results$df_plot.within.2 <- data.frame(actual=pull(df_model2, get(response_var)), prediction=predictions.2)
    results$df_plot.cross.2 <- data.frame(actual=pull(df_model2, get(response_var)), prediction=predictions.2.cross)
    return(results)
}


add_species_zeros <- function(df_list, remove.indv=TRUE) {
    library('stringr')

    all_cols <- c()
    for (df in df_list)
        all_cols <- c(all_cols, colnames(df))

    all_cols <- all_cols[str_detect(all_cols,"s__")]
    all_cols <- unique(all_cols)

    # Add zeros for undetected species
    df_list_mod <- list()
    for (i in seq_along(df_list)) {
        df <- df_list[[i]]
        for (col in all_cols) {
            if (!col%in%colnames(df))
                df[,col] <- 0
        }
        # Put all the columns in the same order
        df <- df[,sort(colnames(df))]
        # Remove the "Individual" column
        if (remove.indv)
            df <- df[,!colnames(df)%in%"Individual"]
        df_list_mod[[i]] <- df
    }
    return(df_list_mod)

}

relative_microbial_age_models <- function(df_model_build, ntree=10000) {
    
    model <- RRF(age ~ ., df_model_build, ntree=ntree)
    predictions <- predict_in_cv(df_model_build, ntree=ntree)
    model_spline_healthy <- smooth.spline(x = df_model_build$age, y = predictions, df=3)
    return(list(model=model, predictions=predictions, model_spline_healthy=model_spline_healthy))
}

relative_microbial_age <- function(df_model_build, df_model_predict, model_components=NULL, ntree=10000) {
    library('RRF')
    if (is.null(model_components)) {
        model_components <- relative_microbial_age_models(df_model_build, ntree)
    }

    if (!identical(df_model_build,df_model_predict))
        predictions_df2 <- predict(model_components$model, df_model_predict)
    else
        predictions_df2 <- model_components$predictions
    df_plot <- data.frame(age=df_model_predict$age, relative_age=predictions_df2 - predict(model_components$model_spline_healthy, df_model_predict$age)$y)
    rownames(df_plot) <- rownames(df_model_predict)
    return(list(df_plot=df_plot, model_components=model_components))
}

### Dissecting the model

calc_average_abundances <- function(df_model, species.oi) {

    vals_healthy <- list()
    for (sp in species.oi) {
        vals <- pull(df_model, get(sp))
        vals_healthy[[sp]] <- mean(vals[vals>0])
    }
    return(vals_healthy)
}

make_blank_model_df <- function(species.use, vals_healthy) {
    df_model_sub <- data.frame(age = 1:2)
    
    for (sp in species.oi) {
        df_model_sub[,sp] <- 0
    }
    
    for (species in species.use) {
        df_model_sub[1:nrow(df_model_sub),species] <- vals_healthy[[species]]
    }

    return(df_model_sub)
}


do_knockout <- function(df_model_sub, species, r2) {
    df_model_sub_copy <- df_model_sub
    df_model_sub_copy[1:nrow(df_model_sub_copy),species] <- 0
    r3 <- relative_microbial_age(df_model_sub, df_model_sub_copy, r2$model_components, ntree=1000)
    return(mean(r2$df_plot$relative_age - r3$df_plot$relative_age))

}


calc_interaction_val <- function(df_model, species.oi, species.test) {
    stopifnot(length(species.test) == 2)
    df_model_sub <- df_model[,c(species.oi,"age")] %>% na.omit

    #Calculate the age model components for future use
    r <- relative_microbial_age(df_model_sub, df_model_sub, ntree=1000)
    model_components <- r$model_components
    
    vals_healthy <- calc_average_abundances(df_model, species.oi)
    
    i <- 1
    max.do <- 100
    df_out <- data.frame()
    for (species.use in combn(species.oi,6,simplify=FALSE)) {
        if (sum(species.test%in%species.use) != 2)
            next
        df_model_dummy <- make_blank_model_df(species.use, vals_healthy)
        r2 <- relative_microbial_age(df_model_sub, df_model_dummy, model_components, ntree=1000)

        # Now do pairwise knockouts
        vals <- c()
        for (species in species.test) {
            vals <- c(vals, do_knockout(df_model_dummy, species, r2))
        }
        vals <- c(vals, do_knockout(df_model_dummy, species.test, r2))

        vals <- data.frame(vals) %>% t  %>% data.frame
        colnames(vals) <- c(species.test, "double")

        df_out <- bind_rows(df_out, vals)
        i <- i + 1
        if (i >= max.do)
            break
    }

    # Calculate the interaction term, based on an additive model
    df_out$diff <- df_out[,1] + df_out[,2] - df_out[,3]

    return(c(mean(df_out$diff),wilcox.test(df_out$diff)$p.value))
}



identify_subspecies <- function(df_model, metric.use="cor", metric.cutoff=0.95, ntree=1000) {
    df_model <- df_model %>% select(-Individual)
    model <- RRF(age ~ ., data=df_model, ntree=ntree)
    df_importance <- get_importance_df(model)
    df_submodels <- build_submodels(df_importance, df_model, n.plot=30, y.var="age", metric ="cor", ntree=ntree)
    metric.use.ratio <- paste0(metric.use,".ratio")

    df_submodels[,metric.use.ratio] <- pull(df_submodels,get(metric.use)) / max(pull(df_submodels,get(metric.use)))
    n.use <- which(pull(df_submodels,get(metric.use.ratio)) >= metric.cutoff)[1]
    species.oi <- head(df_importance, n.use) %>% pull(species)
    return(species.oi)
}


add_missing_species <- function(df_model, species.oi) {
    for (sp in species.oi) {
        if (!sp%in%colnames(df_model))
            df_model[,sp] <- 0
    }
    return(df_model)

}


identify_colonization_points <- function(df_prevalence, species.oi, prevalence.cutoff=0.25) {
    colonization_healthy <- list()
    for (sp in species.oi) {
        df_out_sub <- filter(df_prevalence, (species == sp) & (prevalence >= prevalence.cutoff))
        age_range_use <- df_out_sub %>% pull(age_range) %>% as.character
        vals <- str_split(age_range_use[1],"-")[[1]] %>% as.numeric
        # Midpoint of the range
        colonization_healthy[[sp]] <- sum(vals)/2
    }
    return(colonization_healthy)

}

simplify_name <- function(sp) {
    # Simplifies to Fp
    return(str_extract_all(sp, "_[A-Za-z]{1}")[[1]] %>% str_replace("_","") %>% paste(collapse=""))
}

simplify_name_standard <- function(sp) {
    # Simplifies to F. prausnitzii
    str_replace(sp, "s__([A-Z]{1})[a-z]+","\\1") %>% str_replace("_" ,". ")
}

simplify_name_full<- function(sp) {
    # Simplifies to Faecalibacterium prausnitzii
    str_replace(sp, "s__","") %>% str_replace("_" ," ")
}



check_consecutive <- function(v, n.consecutive.fail=2) {
    for (i in 1:(length(v)-n.consecutive.fail+1)) {
        if (all(v[i:(i+n.consecutive.fail-1)] < 0))
            return(FALSE)
    }
    return(TRUE)
}


summarize_species_detected <- function(species_detected, species.all, n_comb=NULL, total.combs=NULL) {
    if (is.null(total.combs))
        total.combs <- length(combn(species.all,n_comb,simplify=FALSE))
    print(sprintf("total percentage of combinations that pass %f", length(species_detected)/total.combs * 100))

    df_plot <- NULL
    for (sp in species.all) {
        n_found <- lapply(species_detected, function(x) sp%in%x) %>% unlist %>% sum
        frac <- n_found / length(species_detected)
        df_plot <- bind_rows(df_plot, data.frame(sp, frac))
    }
    df_plot <- arrange(df_plot,desc(frac))
    df_plot$sp <- str_replace(df_plot$sp,"s__","") %>% str_replace_all("_"," ")

    df_plot$sp <- factor(df_plot$sp, levels=df_plot$sp, ordered=TRUE)
    return(df_plot)
}
