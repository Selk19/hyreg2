

#' give_id: function to decode which id or observation was classified to which class by the model
#'
#' @description This function can be used to decode the classified classes by the model and see, which id ot observation
#'              was signed to which class
#'
#' @param data a dataframe, which was used to estimate the model
#' @param model a flexmix model object estimtaed using hyreg2 or hyreg2_het
#' @param id one character, indicator for the column holding the ids, this column must be part of the provided dataframe.
#'          the parameter must be specified, if the provided model was estimated under controll for | id
#'
#'
#' @return dataframe of two columns, first column named as provided id column or "observation" if id was not given as
#'          an input. second column named mod_comp indicating the assigned class for this id or observation
#'
#' @details PUT DETAILS HERE
#'
#' @example
#' # estimtae a model using simulated_data_rnorm
#'
#' ### using | id ####
#'formula <- y ~  -1 + x1 + x2 + x3 | id
#'k <- 2
#'stv <- setNames(c(0.2,0.2,0.2,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
#'control <- list(iter.max = 1000, verbose = 4)
#'
#'hyflex_mod <- hyreg2(formula = formula,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type,
#'                     stv = stv,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "cont",
#'                     id_col = "id"
#')

#' # use of function give_id
#' give_id(data = simulated_data_norm,
#' model = hyflex_mod,
#' id = "id")
#'
#'
#'  ### model without control for id during | id ###
#'formula <- y ~  -1 + x1 + x2 + x3
#'k <- 2
#'stv <- setNames(c(0.2,0.2,0.2,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
#'control <- list(iter.max = 1000, verbose = 4)
#'
#'hyflex_mod <- hyreg2(formula = formula,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type,
#'                     stv = stv,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "both"
#')
#' # use of function give_id
#' give_id(data = simulated_data_norm,
#' model = hyflex_mod)
#'
#' @author Svenja Elkenkamp & Kim Rand
#' @export
#'
#'


give_id <- function(data,
                    model,
                    id = NULL) # id must be provided, if model was estimated using | id
  {

    if(!is.list(model)){
      if(!is.null(id)){
        ids_comp <- data.frame(unique(data[,id]),model@cluster)
        colnames(ids_comp) <- c(id,"mod_comp")
      }else{
        data$mod_comp <- model@cluster
        ids_comp <- data.frame(rownames(data),data[,"mod_comp"])
        colnames(ids_comp) <- c("observation","mod_comp")
      }

    }else{
      ids_comp <- model[[3]]
      colnames(ids_comp) <- c(id,"mod_comp")
    }

    return(ids_comp)
}








#' plot_hyreg2: plot function to visualize the classification based on the model
#'
#' @description This function can be used to visualize the classification based on the model for different variables.
#'              ggplot is used.
#'
#' @param data a dataframe, which was used to estimate the model
#' @param x one charachter, column of data to be plotted in x-axis
#' @param y one charachter, column of data to be plotted in y-axis
#' @param id_col_data one charachter, column of data to identify the ids, same as was given in model.
#'            if model was estimated without | id, see details
#' @param id_df_model dataframe of two columns indicating which id belongs to which class,
#'                  first column named id_col_data, second column named mod_comp.
#'                  this variable can easily be filled using the output of give_id().
#' @param colors charachter vector, colors to be used in ggplot, default is NULL - than colors are choosen automalically
#'
#'
#'
#' @return ggplot object visulizing x against y by classes from the model
#'
#' @details
#' id_col_df has to be provided anyway, even if the model was estimated without | id.
#' We recommend to create a new column "observation" in data using the rownames/observationnumbers as charachter values
#' and use this column as input for id_col_data in plot_hyreg2. For id_df_model you can use give_id(data,model,"observation")
#'
#' @example
#' # estimtate a model using simulated_data_rnorm
#'
#' ### using | id ####
#'formula <- y ~  -1 + x1 + x2 + x3 | id
#'k <- 2
#'stv <- setNames(c(0.2,0.2,0.2,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
#'control <- list(iter.max = 1000, verbose = 4)
#'
#'hyflex_mod <- hyreg2(formula = formula,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type,
#'                     stv = stv,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "cont",
#'                     id_col = "id"
#')

#'# plotting the variables x2 against y
#'plot_hyreg(data = simulated_data_norm,
#'           x = "x2",
#'           y = "y",
#'           id_col_data = "id",
#'           id_df_model = give_id(data = simulated_data_norm,
#'                                 model = hyflex_mod,
#'                                 id = "id"))
#'
#'
#'  ### model without control for id during | id ###
#'formula <- y ~  -1 + x1 + x2 + x3
#'k <- 2
#'stv <- setNames(c(0.2,0.2,0.2,1,1),c(colnames(simulated_data_norm)[3:5],c("sigma","theta")))
#'control <- list(iter.max = 1000, verbose = 4)
#'
#'hyflex_mod <- hyreg2(formula = formula,
#'                     data =  simulated_data_norm,
#'                     type =  simulated_data_norm$type,
#'                     stv = stv,
#'                     k = k,
#'                     type_cont = "TTO",
#'                     type_dich = "DCE_A",
#'                     opt_method = "L-BFGS-B",
#'                     control = control,
#'                     latent = "both"
#')
#'
#'# plotting the variables x2 against y
#'
#' simulated_data_norm_plot <- simulated_data_norm
#' simulated_data_norm_plot$observation <- rownames(simulated_data_norm_plot)
#'plot_hyreg(data = simulated_data_norm_plot,
#'           x = "x2",
#'           y = "y",
#'           id_col_data = "observation",
#'           id_df_model = give_id(data = simulated_data_norm_plot,
#'                                 model = hyflex_mod,
#'                                 id = "observation"))
#'
#'
#' @author Svenja Elkenkamp & Kim Rand
#' @importFrom ggplot2 ggplot
#' @export
#'
#'

# include colot option
plot_hyreg <- function(data,
                       x,
                       y,
                       id_col_data,
                       id_df_model,  # you can use give_id() to generate id_df
                       colors = NULL # optional colour vector
){


  colnames(id_df_model) <- c(id_col_data,"mod_comp")
  data <- merge(data, id_df_model, by = id_col_data)

library(ggplot2)
 p <-  ggplot2::ggplot(mapping =aes(x = data[,x], y = data[,y], color = as.character(data$mod_comp)))
 p <- p + geom_point(size = 3, alpha = 0.7)
 p <- p +  labs( title = "classification",
          x = paste0(x),
          y = paste0(y),
          color = "Class")
 p <- p + geom_jitter(width = 0.2, height = 0, size = 2, alpha = 0.7) # verschieben der Punkte leicht zur Seite
 p <- p + theme_minimal()

  # colours
  if (!is.null(colors)) {
    p <- p + scale_color_manual(values = colors)
  }

 return(p)

}





