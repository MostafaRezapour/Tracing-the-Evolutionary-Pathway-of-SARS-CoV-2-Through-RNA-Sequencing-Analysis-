main <- function() {
  tryCatch({
    print("Starting script execution.")
    cat("\014")
    rm(list = ls())
    while (!is.null(dev.list())) dev.off()
    print("Script has started running.")
    library(edgeR)
    library(tcltk)
    library(rstudioapi)
    print("Guideline")
    intro_message <- function() {
      require(tcltk)
      tt <- tktoplevel()
      tkwm.title(tt, "RNA-Seq Analysis Software")
      tcl("wm", "attributes", tt, topmost = TRUE)
      msg <- "To use this software, you need to have the RNA-Seq count data (post-alignment) ready in a CSV file where the rows are genes and columns are samples.
      
    "
      tkgrid(tklabel(tt, text = msg, wraplength = 600, justify = "left"), padx = 20, pady = 20)
      user_choice <- tclVar("continue")
      button_frame <- tkframe(tt)
      continue_btn <- tkbutton(button_frame, text = "Continue", command = function() {
        tclvalue(user_choice) <- "continue"  # Set choice to continue
        tkdestroy(tt)  # Close the window
      })
      exit_btn <- tkbutton(button_frame, text = "Exit and Fix", command = function() {
        tclvalue(user_choice) <- "exit"  # Set choice to exit
        tkdestroy(tt)  # Close the window
      })
      
      # Arrange buttons in the window
      tkgrid(continue_btn, padx = 10, pady = 10)
      tkgrid(exit_btn, padx = 10, pady = 10)
      tkgrid(button_frame)
      tkwait.window(tt)
      if (tclvalue(user_choice) == "exit") {
        cat("Please fix the input CSV file and rerun the script.\n")
        return(FALSE)  
      }

      return(TRUE)  
    }
    
    if (!intro_message()) {
      stop("Execution stopped by user. Fix the input CSV file and rerun the script.")
    }
    
    cat("Script continues...\n")

    load_data <- function(filepath) {
      data <- read.csv(filepath, header = TRUE, check.names = FALSE)
      
      gene_symbols <- data[, 1]  # First column as gene symbols
      counts <- data[, -1]       # All other columns as counts
      
      dge <- DGEList(counts = as.matrix(counts), genes = data.frame(GeneSymbol = gene_symbols))
      o <- order(rowSums(dge$counts), decreasing = TRUE)
      dge <- dge[o, ]
      d <- duplicated(dge$genes$GeneSymbol)
      dge <- dge[!d, ]
      rownames(dge$counts) <- dge$genes$GeneSymbol
      return(dge)
    }
    detect_groups <- function(column_names) {
      groups <- sub(" \\(.*\\)", "", column_names)
      
      group_info <- table(groups)
      return(as.list(group_info))
    }
    
    confirm_groups <- function(group_info) {
      tt <- tktoplevel()
      tkwm.title(tt, "Confirm Groups")
      tkgrid(tklabel(tt, text = "Detected Groups and Sample Counts:"))
      for (group_name in names(group_info)) {
        tkgrid(tklabel(tt, text = sprintf("%s: %d samples", group_name, group_info[[group_name]])))
      }
      user_decision <- tclVar("no")
      tkgrid(tkbutton(tt, text = "Confirm", command = function() {
        tclvalue(user_decision) <- "yes"
        tkdestroy(tt)
      }))
      tkgrid(tkbutton(tt, text = "Decline", command = function() {
        tclvalue(user_decision) <- "no"
        tkdestroy(tt)
      }))
      tkwait.window(tt)
      return(tclvalue(user_decision) == "yes")
    }
    run_analysis <- function() {
      filepath <- file.choose()
      dge <- load_data(filepath)
      column_names <- colnames(dge$counts)
      group_info <- detect_groups(column_names)
      
      if (confirm_groups(group_info)) {
        group_labels <- rep(NA, ncol(dge$counts))
        for (group_name in names(group_info)) {
          group_indices <- which(sub(" \\(.*\\)", "", column_names) == group_name)
          group_labels[group_indices] <- group_name
        }
        dge <- DGEList(counts = dge$counts, group = factor(group_labels))
        keep <- filterByExpr(dge$counts, group = dge$samples$group)
        dge <- dge[keep, , keep.lib.sizes = FALSE]
        print(dge)
        return(dge)
      } else {
        tkmessageBox(title = "Error", message = "Please fix the CSV file and try again.")
        return(NULL)
      }
    }
    dge <- run_analysis()
    table(dge$samples$group)
    dge <- calcNormFactors(dge)
    norm_factors <- dge$samples$norm.factors  # Retrieve normalization factors
    lib_sizes <- dge$samples$lib.size         # Retrieve library sizes
    effective_lib_sizes <- lib_sizes * norm_factors
    tmm_normalized_counts <- cpm(dge, normalized.lib.sizes = TRUE)
    write.csv(tmm_normalized_counts, file = "TMM_normalized_counts.csv", row.names = TRUE)
    select_groups <- function(dge) {
      require(tcltk)  # Ensure the tcltk library is loaded
      group_names <- levels(dge$samples$group)
      tt <- tktoplevel()
      tkwm.title(tt, "Select Baseline and Treatment Groups")
      baseline <- tclVar(group_names[1])
      treatment <- tclVar(if (length(group_names) > 1) group_names[2] else group_names[1])
      tkgrid(tklabel(tt, text = "Select Baseline Group:"))
      baseline_menu <- ttkcombobox(tt, values = group_names, textvariable = baseline, state = "readonly")
      tkgrid(baseline_menu)
      tkgrid(tklabel(tt, text = "Select Treatment Group:"))
      treatment_menu <- ttkcombobox(tt, values = group_names, textvariable = treatment, state = "readonly")
      tkgrid(treatment_menu)
      confirm_btn <- tkbutton(tt, text = "Confirm", command = function() {
        if (tclvalue(baseline) == tclvalue(treatment)) {
          tkmessageBox(title = "Error", message = "Baseline and Treatment groups must be different.")
        } else {
          tkdestroy(tt)
        }
      })
      tkgrid(confirm_btn)
      tkwait.window(tt)
      baseline_group_size <- sum(dge$samples$group == tclvalue(baseline))
      treatment_group_size <- sum(dge$samples$group == tclvalue(treatment))
      print(paste("Size of Baseline group (", tclvalue(baseline), "): ", baseline_group_size, sep = ""))
      print(paste("Size of Treated group (", tclvalue(treatment), "): ", treatment_group_size, sep = ""))
      return(list(baseline = tclvalue(baseline), treated = tclvalue(treatment)))
    }
    perform_glm_qlf <- function(dge, selected_groups) {
      dge_subset <- dge[, dge$samples$group %in% c(selected_groups$baseline, selected_groups$treated)]
      dge_subset$samples$group <- relevel(dge_subset$samples$group, ref = selected_groups$baseline)
      design <- model.matrix(~ group, data = dge_subset$samples)
      dge_subset <- estimateDisp(dge_subset, design)
      fit <- glmQLFit(dge_subset, design)
      qlf <- glmQLFTest(fit, coef = 2)  # coef = 2 compares the treated group against the baseline group
      top_results <- topTags(qlf, n = Inf)
      print(head(top_results$table))
      write.csv(top_results$table, file = "differential_expression_results.csv", row.names = TRUE)
      plotMD(qlf, main = "MA Plot", xlab = "Average Log CPM", ylab = "Log Fold Change")
      abline(h = c(-1, 1), col = "blue")
      return(top_results)
    }
    ask_comparisons <- function() {
      require(tcltk)
      tt <- tktoplevel()
      tkwm.title(tt, "Number of Comparisons")
      tkgrid(tklabel(tt, text = "How many comparisons would you like to perform?"))
      num_comparisons <- tclVar("1")  # Default value
      entry <- tkentry(tt, textvariable = num_comparisons, width = 10)
      tkgrid(entry)
      confirm_btn <- tkbutton(tt, text = "Confirm", command = function() {
        tkdestroy(tt)
      })
      tkgrid(confirm_btn)
      tkwait.window(tt)
      return(as.integer(tclvalue(num_comparisons)))
    }
    num_comparisons <- ask_comparisons()
    for (i in 1:num_comparisons) {
      cat("\nPerforming comparison", i, "of", num_comparisons, "\n")
      selected_groups <- select_groups(dge)  # Select baseline and treated groups
      results <- perform_glm_qlf(dge, selected_groups)
      filename <- paste0("DE_", selected_groups$baseline, "_vs_", selected_groups$treated, ".csv")
      results_table <- results$table
      results_table <- data.frame(GeneSymbol = rownames(results_table), results_table)
      write.csv(results_table, file = filename, row.names = FALSE)
      cat("Differential expression results saved to:", filename, "\n")
    }
    print("Script completed successfully.")
  }, error = function(e) {
    print(sprintf("An error occurred: %s", e$message))
  })
}
# Call the main function to execute
main()
