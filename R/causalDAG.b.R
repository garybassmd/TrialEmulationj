#' @importFrom jmvcore .
causalDAGClass <- R6::R6Class(
    "causalDAGClass",
    inherit = causalDAGBase,
    private = list(
        .run = function() {
            if (is.null(self$options$exposure) || is.null(self$options$outcome) || !nzchar(trimws(self$options$dagText)))
                return()
            fitted <- tryCatch(private$.analyseDAG(), error = function(e) e)
            if (inherits(fitted, "error"))
                jmvcore::reject(paste0("The DAG could not be validated: ", conditionMessage(fitted)))
            private$.populate(fitted)
        },

        .analyseDAG = function() {
            if (!requireNamespace("dagitty", quietly = TRUE))
                stop("The dagitty dependency is unavailable.")
            text <- trimws(self$options$dagText)
            if (!grepl("^(dag|mag|pdag|pag)\\s*\\{", text, ignore.case = TRUE))
                text <- paste0("dag { ", text, " }")
            graph <- dagitty::dagitty(text, layout = TRUE)
            if (!dagitty::isAcyclic(graph))
                stop("The graph contains a directed cycle. Causal DAGs must be acyclic.")
            exposure <- self$options$exposure
            outcome <- self$options$outcome
            nodes <- names(graph)
            if (!exposure %in% nodes)
                stop("Exposure node '", exposure, "' is not present in the graph. Available nodes: ", paste(nodes, collapse = ", "), ".")
            if (!outcome %in% nodes)
                stop("Outcome node '", outcome, "' is not present in the graph. Available nodes: ", paste(nodes, collapse = ", "), ".")
            if (identical(exposure, outcome))
                stop("Exposure and outcome must be different nodes.")
            dagitty::exposures(graph) <- exposure
            dagitty::outcomes(graph) <- outcome

            sets <- dagitty::adjustmentSets(
                graph,
                exposure = exposure,
                outcome = outcome,
                type = "minimal",
                effect = self$options$effect,
                max.results = self$options$maxSets
            )
            selected <- unique(self$options$selectedAdjusters)
            selected <- selected[!is.na(selected) & nzchar(selected)]
            selectedInGraph <- intersect(selected, nodes)
            selectedValid <- if (identical(self$options$effect, "total")) {
                tryCatch(
                    dagitty::isAdjustmentSet(graph, selectedInGraph, exposure, outcome),
                    error = function(e) FALSE
                ) && length(selectedInGraph) == length(selected)
            } else {
                length(selectedInGraph) == length(selected) && any(vapply(
                    sets,
                    function(set) setequal(selectedInGraph, unname(unlist(set, use.names = FALSE))),
                    logical(1)
                ))
            }

            edges <- dagitty::edges(graph)
            directed <- edges[as.character(edges$e) == "->", , drop = FALSE]
            incoming <- table(directed$w)
            colliders <- names(incoming)[incoming >= 2]
            descendants <- dagitty::descendants(graph, exposure, proper = TRUE)
            minimalNodes <- unique(unlist(sets, use.names = FALSE))
            selectedRows <- lapply(selected, function(variable) {
                if (!variable %in% nodes) {
                    status <- "Not in DAG"
                    explanation <- "Add the node and its causal arrows or remove it from the adjustment model."
                } else if (variable %in% descendants && identical(self$options$effect, "total")) {
                    status <- "Do not adjust"
                    explanation <- "This node is downstream of the exposure and is not eligible for total-effect adjustment."
                } else if (variable %in% colliders) {
                    status <- "Collider warning"
                    explanation <- "This node has two or more incoming arrows; conditioning may open a noncausal path."
                } else if (variable %in% minimalNodes) {
                    status <- "Supported"
                    explanation <- "This node appears in at least one minimally sufficient adjustment set."
                } else {
                    status <- "Review"
                    explanation <- "This node is not required by the displayed minimal sets. Confirm its role and avoid unnecessary adjustment."
                }
                data.frame(variable = variable, status = status, explanation = explanation, stringsAsFactors = FALSE)
            })
            selectedTable <- if (length(selectedRows) == 0L) data.frame(variable = character(), status = character(), explanation = character()) else do.call(rbind, selectedRows)

            coordinates <- dagitty::coordinates(graph)
            if (is.null(coordinates$x) || any(!is.finite(coordinates$x))) {
                graph <- dagitty::graphLayout(graph)
                coordinates <- dagitty::coordinates(graph)
            }
            nodeData <- data.frame(
                node = names(coordinates$x),
                x = as.numeric(coordinates$x),
                y = as.numeric(coordinates$y),
                stringsAsFactors = FALSE
            )
            nodeData$role <- "Other"
            nodeData$role[nodeData$node %in% dagitty::latents(graph)] <- "Unobserved"
            nodeData$role[nodeData$node %in% selected] <- "Selected adjuster"
            nodeData$role[nodeData$node == exposure] <- "Exposure"
            nodeData$role[nodeData$node == outcome] <- "Outcome"
            edgeData <- edges[, c("v", "w", "e"), drop = FALSE]
            edgeData <- merge(edgeData, nodeData[, c("node", "x", "y")], by.x = "v", by.y = "node", all.x = TRUE)
            names(edgeData)[names(edgeData) %in% c("x", "y")] <- c("xStart", "yStart")
            edgeData <- merge(edgeData, nodeData[, c("node", "x", "y")], by.x = "w", by.y = "node", all.x = TRUE)
            names(edgeData)[names(edgeData) %in% c("x", "y")] <- c("xEnd", "yEnd")

            independencies <- tryCatch(dagitty::impliedConditionalIndependencies(graph), error = function(e) list())
            fingerprint <- .te_specification_fingerprint(list(
                dag = text, exposure = exposure, outcome = outcome,
                effect = self$options$effect, selected = selected
            ))
            list(
                graph = graph,
                nodes = nodeData,
                edges = edgeData,
                sets = sets,
                selected = selectedTable,
                selectedValid = selectedValid,
                independencies = independencies,
                exposure = exposure,
                outcome = outcome,
                colliders = colliders,
                descendants = descendants,
                fingerprint = fingerprint
            )
        },

        .populate = function(fitted) {
            rows <- list(
                c("Graph type", dagitty::graphType(fitted$graph)),
                c("Nodes", nrow(fitted$nodes)),
                c("Edges", nrow(fitted$edges)),
                c("Exposure", fitted$exposure),
                c("Outcome", fitted$outcome),
                c("Effect", if (self$options$effect == "total") "Total" else "Direct"),
                c("Minimal adjustment sets", length(fitted$sets)),
                c("Selected set valid", if (fitted$selectedValid) "Yes" else "No / incomplete")
            )
            for (i in seq_along(rows))
                self$results$summary$setRow(rowNo = i, values = list(measure = rows[[i]][1], value = rows[[i]][2]))
            if (length(fitted$sets) == 0L) {
                self$results$adjustmentSets$addRow(rowKey = "set1", values = list(number = 1, variables = "No adjustment required, or effect not identifiable by covariate adjustment"))
            } else {
                for (i in seq_along(fitted$sets)) {
                    variables <- unname(unlist(fitted$sets[[i]], use.names = FALSE))
                    self$results$adjustmentSets$addRow(rowKey = paste0("set", i), values = list(
                        number = i,
                        variables = if (length(variables) == 0L) "Empty set (no adjustment)" else paste(variables, collapse = ", ")
                    ))
                }
            }
            if (nrow(fitted$selected) > 0L) {
                for (i in seq_len(nrow(fitted$selected)))
                    self$results$selected$addRow(rowKey = paste0("selected", i), values = as.list(fitted$selected[i, ]))
            }
            independenceText <- if (length(fitted$independencies) == 0L) "No nontrivial conditional independencies were returned." else paste(capture.output(print(fitted$independencies)), collapse = "\n")
            self$results$independencies$setContent(independenceText)
            self$results$dagPlot$setState(list(nodes = fitted$nodes, edges = fitted$edges))
            self$results$notes$setContent(paste0("- ", c(
                "The graph is a prespecified causal model, not a structure learned from this dataset.",
                "A valid adjustment set follows only if the graph is causally correct and all nodes designated as observed are measured adequately.",
                "Do not choose covariates because they predict the outcome or improve balance alone. Use temporal and causal knowledge before inspecting treatment effects.",
                if (length(fitted$colliders) > 0L) paste0("Potential collider nodes based on incoming arrows: ", paste(fitted$colliders, collapse = ", "), ".") else "No nodes with two or more incoming directed arrows were detected.",
                if (length(fitted$descendants) > 0L) paste0("Exposure descendants: ", paste(fitted$descendants, collapse = ", "), ". These are generally excluded from total-effect adjustment.") else "No exposure descendants were detected.",
                "DAGs do not by themselves verify consistency, positivity, measurement quality, or absence of selection bias.",
                if (self$options$effect == "direct") "For direct effects, the selected-set check is conservative and requires an exact match to one displayed minimal direct-effect set." else NULL,
                paste0("Specification fingerprint: ", substr(fitted$fingerprint, 1, 16))
            ), collapse = "\n"))
            self$results$syntax$setContent(paste0(
                as.character(fitted$graph),
                "\n\n# Minimal adjustment sets\n",
                paste(capture.output(print(fitted$sets)), collapse = "\n"),
                "\n# Specification SHA-256: ", fitted$fingerprint
            ))
        },

        .dagPlot = function(image, ggtheme, theme, ...) {
            if (is.null(image$state))
                return(FALSE)
            nodes <- image$state$nodes
            edges <- image$state$edges
            directed <- edges[as.character(edges$e) == "->", , drop = FALSE]
            other <- edges[as.character(edges$e) != "->", , drop = FALSE]
            plot <- ggplot2::ggplot()
            if (nrow(other) > 0L)
                plot <- plot + ggplot2::geom_segment(data = other, ggplot2::aes(x = xStart, y = yStart, xend = xEnd, yend = yEnd), linetype = "dashed", colour = "#666666", linewidth = .7)
            if (nrow(directed) > 0L)
                plot <- plot + ggplot2::geom_segment(data = directed, ggplot2::aes(x = xStart, y = yStart, xend = xEnd, yend = yEnd), arrow = grid::arrow(length = grid::unit(0.16, "inches"), type = "closed"), lineend = "round", colour = "#555555", linewidth = .7)
            plot +
                ggplot2::geom_point(data = nodes, ggplot2::aes(x = x, y = y, fill = role), shape = 21, size = 8, colour = "#333333", stroke = .8) +
                ggplot2::geom_text(data = nodes, ggplot2::aes(x = x, y = y, label = node), size = 3.5) +
                ggplot2::scale_fill_manual(values = c("Exposure" = "#E69F00", "Outcome" = "#0072B2", "Selected adjuster" = "#009E73", "Unobserved" = "#999999", "Other" = "#F2F2F2")) +
                ggplot2::coord_equal() +
                ggplot2::labs(x = NULL, y = NULL, fill = "Node role") +
                ggtheme +
                ggplot2::theme(axis.text = ggplot2::element_blank(), axis.ticks = ggplot2::element_blank(), panel.grid = ggplot2::element_blank())
        }
    )
)
