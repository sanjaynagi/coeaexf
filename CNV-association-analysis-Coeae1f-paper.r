library(data.table)
library(magrittr)
library(stringr)
library(glmmTMB)

git.folder <- 'https://raw.githubusercontent.com/vigg-lstm/GAARD_work/v2.0'
get.git.path <- function(filepath){
	filepath <- paste(filepath, collapse = '/')
	paste(git.folder, filepath, sep = '/')
}

source(get.git.path('shared_functions/R_plotting.r'))

study.ids.table <- fread(get.git.path('data/study_ids.csv'))
study.ids.table[, study_id := paste('v3.2_', study_id, sep = '')]
meta <- fread(get.git.path('data/combined/all_samples.samples.meta.csv'), key = 'sample_id')
meta.sample.names <- meta$partner_sample_id
phen <- fread(get.git.path('data/combined/sample_phenotypes.csv'), key = 'specimen')[sort(meta.sample.names), ]

# Load the sib groups information
sib.groups <- fread(get.git.path('NGSrelate/full_relatedness/sib_group_table.csv'), sep = '\t')
sibs.to.remove <- sib.groups[keep == F, sample.name]

# Identify the males and add them to the list of samples to remove
males <- meta[sex_call == 'M', partner_sample_id]
samples.to.remove <- c(sibs.to.remove, males)

# Now get the modal CNVs by gene
load.modal.cnv.table <- function(study.id){
	filepath <- get.git.path(c('CNV_analysis/Ag1000G_CNV_data', 
	                           study.id, 
	                           'modal_CNVs/modal_copy_number_gambcolu.csv'))
	modal.cnv.table <- fread(filepath)
	these.sample.names <- meta[modal.cnv.table$V1, partner_sample_id]
	modal.cnv.table$V1 <- these.sample.names
	colnames(modal.cnv.table)[1] <- 'sample.id'
	modal.cnv.table
}

modal.copy.number <- lapply(unique(study.ids.table$study_id), load.modal.cnv.table) %>%
                     rbindlist() %>%
                     .[!(sample.id %in% samples.to.remove)] %>%
                     .[high_variance == F, -c('sex_call', 'high_variance')]

modal.CNV.table <- copy(modal.copy.number) %>%
                   .[, colnames(.)[-1] := lapply(.SD, `>`, 0), .SDcols = colnames(.)[-1]]

# For a smaller table, we only keep genes where a CNV was present. 
no.cnv <- names(which(apply(modal.CNV.table[, -1], 2, sum, na.rm = T) < 1))
modal.CNV.table <- modal.CNV.table[, -c(..no.cnv)]

sample.names <- modal.copy.number$sample.id

# Get the size of each population-year combination
population.popsize <- interaction(phen[sample.names, .(location, species)]) %>%
                      droplevels %>%
                      table

population.modal.CNV.counts <- aggregate(modal.CNV.table[, -c('sample.id')],
                                         phen[modal.CNV.table$sample.id, .(location, species)], 
                                         function(x) sum(x)
)

# Calculate the frequency of CNVs in each population by year
population.modal.CNVs <- aggregate(modal.CNV.table[, -c('sample.id')],
                                   phen[modal.CNV.table$sample.id, .(location, species)], 
                                   function(x) sum(x) / length(x)
)

genes.of.interest <- c('Coeae1f', 'Ace1')

# Get a table associating gene names with IDs
gene.table <- fread(get.git.path('CNV_analysis/Ag1000G_CNV_data/gene_annotation_fullgenetable.csv'), 
                    key = 'Gene_stable_ID', check.names = T)
# Coeae1f is absent from the gff, so add it manually
gene.table['AGAP006227', Gene.name := 'COEAE1F']
gene.name.conversion <- gene.table[Gene.name %in% toupper(genes.of.interest), 
                                   .(Gene.id = Gene_stable_ID, Gene.name = str_to_title(Gene.name))]

# We shrink the modal copy number table to be just the genes we are interested in. 
modal.copy.number <- modal.copy.number[, c('sample.id', gene.name.conversion$Gene.id), with = F] %>%
                     merge(phen[sample.names],
                           by.x = 'sample.id', by.y = 'specimen',
                           all = F
                     ) %>%
                     .[, phenotype := as.factor(phenotype)] %>%
                     setnames(gene.name.conversion$Gene.id, gene.name.conversion$Gene.name) %>%
                     setkey(location, insecticide)


# A function that will round to x significant digits, unless the number is small, in which case it 
# just takes 1 significant digit.
signif2 <- function(x, digits, threshold = 0.1){
	small.x <- x < threshold 
	x[small.x] <- signif(x[small.x], 1)
	x[!small.x] <- signif(x[!small.x], digits)
	x
}

contable <- function(genes, shrink = T, include.plot = T, add.sample.size = F, 
                     star.small.samples = T, use.agaps = F, 
                     species.colours = c(coluzzii = 'coral3', gambiae = 'dodgerblue3'), 
                     text.cell.cex = 0.65, pop.cex = 0.5, gene.cex = 0.5, ...){
	
	cnvs.by.population <- population.modal.CNVs
	#
	if (genes[1] %in% gene.table$Gene_stable_ID)
		gene.names <- gene.table[genes, Gene.name]
	else {
		gene.names <- character()
		for (i in 1:length(genes)){
			gene <- genes[i]
			if (toupper(gene) %in% toupper(gene.table$Gene.name)){
				gene.name <- gene.table$Gene.name[toupper(gene.table$Gene.name) == toupper(gene)]
				if (length(gene.name) != 1) 
					stop('There should be only one gene matching the gene name.')
				gene.names[i] <- str_to_title(gene.name)
				genes[i] <- gene.table[Gene.name == gene.name, Gene_stable_ID]
			}
			else
				stop(paste('Could not find the requested gene (', gene, ') in gene.table.', sep = ''))
		}
	}
	population <- droplevels(interaction(cnvs.by.population[, c('location', 'species')]))
	output.table.rownames <- levels(population)
	output.table.colnames <- if(use.agaps) genes else gene.names
	output.table <- matrix(NA, length(output.table.rownames), length(output.table.colnames), dimnames = list(output.table.rownames, output.table.colnames))
	gene.has.cnv <- genes %in% colnames(cnvs.by.population)
	for (ro in output.table.rownames){
		which.row <- which(population == ro)
		if (length(which.row) != 1)
			output.table[ro, ] <- NA
		else {
			for (i in 1:ncol(output.table)){
				co <- output.table.colnames[i]
				if (gene.has.cnv[i])
					output.table[ro, co] <- as.numeric(cnvs.by.population[which.row, genes[i]])
				else
					output.table[ro, co] <- 0
			}
		}
	}
	if (shrink){
		# We keep the rows where at least one entry is > 0. We then remove the columns where all the remaining
		# rows are NA.
		which.rows <- apply(output.table, 1, sum, na.rm = T) > 0
		output.table <- output.table[which.rows, , drop = F]
		which.columns <- apply(output.table, 2, sum, na.rm = T) > 0
		output.table <- output.table[, which.columns, drop = F]
	}
	if (include.plot){
		if (is.null(names(species.colours)) | !any(names(species.colours) %in% c('gambiae', 'coluzzii')))
			names(species.colours) <- c('gambiae', 'coluzzii')
		par(mar = c(0,3,3,0), ...)
		plot(c(0,ncol(output.table)), -c(0,nrow(output.table)), type = 'n', bty = 'n', xaxt = 'n', yaxt = 'n', xlab = '', ylab = '')
		x.matrix <- matrix(1:ncol(output.table), nrow(output.table), ncol(output.table), byrow = T)
		y.matrix <- matrix(1:nrow(output.table), nrow(output.table), ncol(output.table))
		col.matrix <- matrix(rgb(0.8, 0.8, 0.8), nrow(output.table), ncol(output.table))
		popsize <- population.popsize[output.table.rownames]
		text.col.matrix <- c('black', 'white')[(output.table > 0.5) + 1]
		if (add.sample.size){
			text.cell.matrix <- sub('.*\\.', '.', paste(signif2(output.table, 2), ' (', popsize, ')', sep = ''))
			text.cell.matrix[grepl('NA', text.cell.matrix)] <- NA
		}
		else if (star.small.samples){
			text.cell.matrix <- sub('.*\\.', '.', paste(signif2(output.table, 2), c('', ' *', ' **')[(popsize < 5) + (popsize == 1) + 1], sep = ''))
			text.cell.matrix[grepl('NA', text.cell.matrix)] <- NA
		}
		else {
			text.cell.matrix <- signif2(output.table, 2)
		}
		for (i in 1:nrow(output.table)){
			this.row <- output.table[i, ]
			if (any(!is.na(this.row))){
				this.species <- sub('^.*\\.', '', rownames(output.table)[i])
				col.matrix[i, !is.na(this.row)] <- lighten.col(species.colours[this.species], sqrt(this.row[!is.na(this.row)]))
			}
		}
		rect(x.matrix - 1, 1 - y.matrix, x.matrix, - y.matrix, col = col.matrix, border = 'white', lwd = text.cell.cex)
		text(x.matrix - 0.5, 0.5 - y.matrix, text.cell.matrix, cex = text.cell.cex, col = text.col.matrix, font = 2)
		text(x.matrix[1,] - 0.5, (par('usr')[4] - par('usr')[3])/100, srt = 45, adj = c(0,0), colnames(output.table), xpd = NA, cex = gene.cex)
		country.names <- sub('\\.[^.]*$', '', rownames(output.table))
		species <- sub('.*\\.', '', rownames(output.table))
		rect(x.matrix[1,1]-1, 1 - y.matrix, par('usr')[1] - (par('usr')[2] - par('usr')[1])*2, -y.matrix, col = sapply(species.colours[species], lighten.col, 0.4), border = NA, xpd = NA)
		text(x.matrix[1,1]-1, 0.5 - y.matrix[,1], paste(country.names, '  \n', species, '  ', sep = ''), adj = 1, cex = pop.cex, xpd = NA)
	}
	invisible(output.table)
}

png('Modal_CNVs.png', width = 0.7, height = 1.3, units = 'in', res = 300)
par(family = 'Arial')
CNV.frequencies <- contable(grep('Coeae', genes.of.interest, value = T), 
                            text.cell.cex = 0.35,
                            pop.cex = 0.35,
                            gene.cex = 0.35,
                            mai = c(0,0.3,0.18,0.02),
                            shrink = F
)
dev.off()

glm.up <- function(input.table, list.of.markers = markers, rescolumn = 'phenotype', control.for = character(), glm.function = NULL, verbose = T){
	# Check whether the markers and random effects are present in the data.frame
	if (sum(list.of.markers %in% colnames(input.table)) != length(list.of.markers))
		stop('Some of the requested markers were not found in genotypes table.')
	if (sum(control.for %in% colnames(input.table)) != length(control.for))
		stop('At least one random effect was not found in genotypes table.')
	if (!(rescolumn %in% colnames(input.table)))
		stop('Resistance column not found in genotypes table.')
	# Remove any requested random effects that have only 1 level
	random.effects <- character()
	for (this.control in control.for){
		level.counts <- tapply(input.table[,this.control], input.table[,this.control], length)
		if (max(level.counts, na.rm = T) < sum(level.counts, na.rm = T)) 
			random.effects <- c(random.effects, this.control)
		else
			cat('Removing random effect ', this.control, ' was removed because it is invariable.')
	}
	# If you did not set a glm.function, decide which glm function you are going to use, based on whether mixed 
	# modelling will be necessary
	if (is.null(glm.function))
		glm.function <- ifelse(length(random.effects) > 0, 'glmer', 'glm')
	if (glm.function == 'glm'){
		# Create the random effect string, which will be empty if we have no random effects
		random.effects.string <- ifelse(length(random.effects) > 0, paste(' +', paste(random.effects, collapse = ' + ', sep = '')), '')
		# Set the name of the column containing the P.value in the anova function (which is different depending
		# on the glm function you use
		P.val.column <- 'Pr(>Chi)'
	}
	else if (glm.function %in% c('glmer', 'glmmTMB')){
		random.effects.string <- ifelse(length(random.effects) > 0, paste(' +', paste('(1|', random.effects, ')', collapse = ' + ', sep = '')), '')
		P.val.column <- 'Pr(>Chisq)'
	}
	if (verbose)
		cat('Using the following string to control for confounding factors: ', random.effects.string, '\n', sep = '')
	# We remove markers for which there is no variation in the dataset or for which some alleles are too rare. 
	if (verbose)
		cat('\nDetecting invariable and nearly invariable markers.\n')
	kept.markers <- character()
	invariable.markers <- character()
	for (this.marker in list.of.markers){
		allele.counts <- tapply(input.table[,this.marker], input.table[,this.marker], length)
		if (max(allele.counts, na.rm = T) <= (sum(allele.counts, na.rm = T) - 2)) 
			kept.markers <- c(kept.markers, this.marker)
		else
			invariable.markers <- c(invariable.markers, this.marker)
	}
	if (length(kept.markers) == 0)
		stop('Fail. None of the markers provided were sufficiently variable.')
	# We check whether there are any ordered factors and recode them as numeric
	converted.table <- input.table
	has.ordered.factors <- F
	for (this.marker in kept.markers){
		if ('ordered' %in% class(converted.table[, this.marker])){
			if (verbose)
				cat('Converting ordered factor ', this.marker, ' to numeric.\n', sep = '')
			converted.table[, this.marker] <- as.numeric(converted.table[, this.marker])
			has.ordered.factors <- T
		}
	}
	# We do the glm analysis directly on the table from the global environment rather than the argument, this 
	# way the table that was used is recorded in the output. If we had to convert the ordered factors, then
	# we are forced to use a new table
	if (has.ordered.factors){
		working.table.name <- make.names(paste(deparse(substitute(input.table)), '_numeric_conversion', sep = ''))
		eval(parse(text = paste(working.table.name, '<- converted.table')))
	}
	else
		working.table.name <- deparse(substitute(input.table))
	
	# For each marker, we calculate its pseudo-R2 and P-value compared to the null model.
	if (verbose)
		cat('\nAnalysing markers independently.\n')
	individual.markers <- data.frame(P = numeric(), pseudo.R2 = numeric())
	for (this.marker in kept.markers){
		# Remove the Na values for this marker
		this.table <- converted.table[!is.na(converted.table[,this.marker]),]
		# Build the model 
		this.model <- eval(parse(text = paste(glm.function, '(', rescolumn, ' ~ ', this.marker, random.effects.string, ', data = this.table, family = binomial)', sep = '')))
		# Build the null model
		this.null.model <- eval(parse(text = paste(glm.function, '(', rescolumn, ' ~ 1 ', random.effects.string, ', data = this.table, family = binomial)', sep = '')))
		# Get the stats
		this.p <- anova(this.model, this.null.model, test = 'Chisq')[[P.val.column]][2]
		# Report pseudo Rsquared if we used GLM and if we have the modEvA package
		if (('modEvA' %in% (.packages())) & (glm.function == 'glm'))
			this.pseudo.r2 <- mean(unlist(RsqGLM(this.model)))
		else
			this.pseudo.r2 <- NA
		individual.markers[this.marker, ] <- c(this.p, this.pseudo.r2)
	}
	
	# We now build the null model and add markers one by one until all markers are significant
	working.markers <- character()
	# We'll keep track of the markers with perfect correlation
	correlated.markers <- character()
	if (verbose)
		cat('\nRunning commentary on model optimisation:\n\n')
	while(length(working.markers) < length(kept.markers)){
		# Build the model using the working.markers
		if (length(working.markers)){
			old.model.text <- paste(glm.function, '(', rescolumn, ' ~ ', paste(working.markers, collapse = ' + '), random.effects.string, ', data = ', working.table.name, ', family = binomial)', sep = '')
			old.model <- eval(parse(text = old.model.text))
			# Check the remaining markers as some of them may become monomorphic when NAs from the current marker are 
			# taken into account
			for (this.marker in setdiff(kept.markers, working.markers)){
				markers.subset.genotypes <- input.table[, c(working.markers, this.marker), drop = F]
				markers.subset.genotypes <- markers.subset.genotypes[complete.cases(markers.subset.genotypes), , drop = F]
				number.of.alleles <- unique(markers.subset.genotypes[ ,this.marker])
				if (length(number.of.alleles) < 2) {
					if (verbose){
						cat('Removing marker ', this.marker, ' as it has no variation left when NAs at previously added ',
							'loci are removed.\n\n', sep = '')
					}
					kept.markers <- setdiff(kept.markers, this.marker)
				}
			}
			# If we have removed all the remaining markers, we quit the loop and report the final model. 
			if (length(working.markers) == length(kept.markers)){
				final.model <- old.model
				if (verbose){
					cat('\tNo further markers are significant, keeping final model:\n')
					print(final.model)
				}
				break
			}
		}
		else{
			old.model.text <- paste(glm.function, '(', rescolumn, ' ~ 1 ', random.effects.string, ', data = ', working.table.name, ', family = binomial)', sep = '')
			old.model <- eval(parse(text = old.model.text))
		}
		if (verbose)
			cat('Building model:', old.model.text, '\n')
		p.values <- numeric()
		for (this.marker in setdiff(kept.markers, working.markers)){
			new.model <- eval(parse(text = paste('update(old.model, .~.+', this.marker, ')', sep = '')))
			# Check that the new model doesn't have fewer rows than the old one (if the added marker had unique NAs)
			if (length(fitted(new.model)) == length(fitted(old.model))){
				this.p.value <- anova(old.model, new.model, test = 'Chisq')[[P.val.column]][2]
			}
			# Otherwise, we need to rebuild the models with those samples removed.
			else{
				if (has.ordered.factors)
					reduced.input.table <- converted.table[!is.na(converted.table[,this.marker]),]
				else 
					reduced.input.table <- input.table[!is.na(input.table[,this.marker]),]
				if (length(working.markers))
					temp.old.model <- eval(parse(text = paste(glm.function, '(', rescolumn, ' ~ ', paste(working.markers, collapse = ' + '), random.effects.string, ', data = reduced.input.table, family = binomial)', sep = '')))
				else 
					temp.old.model <- eval(parse(text = paste(glm.function, '(', rescolumn, ' ~ 1 ', random.effects.string, ', data = reduced.input.table, family = binomial)', sep = '')))
				new.model <- eval(parse(text = paste('update(temp.old.model, .~.+', this.marker, ')', sep = '')))
				this.p.value <- anova(temp.old.model, new.model, test = 'Chisq')[[P.val.column]][2]
			}
			# If the p.value was NA, then there is perfect correlation or something else was wrong. Set the p-value 
			# to Inf and move on to the next marker
			if (is.na(this.p.value)){
				if (verbose)
					cat('\tCould not calculate p-value for marker ', this.marker, '.\n\n', sep = '')
				p.values[this.marker] <- Inf
			}
			else{
				p.values[this.marker] <- this.p.value
			}
		}
		if (min(p.values) <= 0.05){
			# Add the lowest significant p-value
			if (verbose)
				cat('\tAdding marker ', names(p.values)[which.min(p.values)], ' as the lowest significant marker (P = ', min(p.values), ').\n\n', sep = '')
			marker.to.add <- names(p.values)[which.min(p.values)]
			working.markers <- c(working.markers, marker.to.add)
			if (length(working.markers) == length(kept.markers)){
				# If all markers have been added, then we have the final model
				final.model <- new.model
				if (verbose)
					cat('\tNo markers left to add.\n')
			}
		}
		else {
			final.model <- old.model
			if (verbose){
				cat('\tNo further markers are significant, keeping final model:\n')
				print(final.model)
			}
			break
		}
	}
	if (verbose){
		if (length(working.markers) == 0)
			cat('Final model was the null model.\n\n')
		else 
			cat('Final model contained ', length(working.markers), ' parameters: ', paste(working.markers, collapse = ','), '.\n\n', sep = '')
	}
	# Now get the p-values and pseudo R-squared value for all the variables in the final model, when added as the 
	# last variable
	if (length(working.markers) > 0){
		deviance.effect <- numeric()
		final.p.values <- numeric()
		for (this.marker in working.markers){
			# Remove the Na values for this marker
			this.table <- converted.table[!is.na(converted.table[,this.marker]),]
			# Build the model 
			reduced.final.model <- update(final.model, data = this.table)
			reduced.model <- eval(parse(text = paste('update(reduced.final.model, .~.-', this.marker, ')', sep = '')))
			# Get the stats
			dev.eff <- deviance(reduced.model) - deviance(final.model)
			deviance.effect[this.marker] <- ifelse(is.null(dev.eff), NA, dev.eff)
			final.p.values[this.marker] <- anova(reduced.model, reduced.final.model, test = 'Chisq')[[P.val.column]][2]
		}
	}
	else {
		final.p.values <- NA
		deviance.effect <- NA
	}
	cat('\n')
	#
	final.model.sig <- data.frame(P = final.p.values, deviance = deviance.effect)
	print(final.model.sig)
	list(invariable.markers = invariable.markers, 
	     correlated.markers = correlated.markers,
	     sig.alone = individual.markers,
	     final.model = final.model,
	     final.sig = final.model.sig)
}

# Let's test things using copy number in each population independently
cat('\n\t#######################')
cat('\n\t# Copy number testing #')
cat('\n\t#######################\n')

genes.to.model <- c('Coeae1f','Ace1')
# CNVs in Coeae1f are only found in Baguida and Obuasi
association.tests <- list()
for (insecticide in c('Delta', 'PM')) {
	association.tests[[insecticide]] <- list()
	for (location in c('Baguida', 'Obuasi')){
		# For some reason, this doesn't work if I do the indexing directly. I have to create this object first
		index <- list(location, insecticide)
		cat('\n\tModal copy number test for ', paste(index, collapse = '_'), ':\n\n', sep = '')
		test.table <- as.data.frame(modal.copy.number[index])
		test.result <- glm.up(test.table, genes.to.model, 'phenotype')
		print(test.result)
		association.tests[[insecticide]][[location]] <- test.result
	}
}

# Let's try a glm that includes population as a random factor. 
delta.table <- as.data.frame(modal.copy.number[(insecticide == 'Delta') & (location %in% c('Baguida', 'Obuasi'))])
cat('\n\nAll populations Delta:\n')
association.tests[['Delta']][['combined']] <- glm.up(delta.table, genes.to.model, 'phenotype', control.for = 'location', glm.function = 'glmmTMB')
print(association.tests[['Delta']][['combined']])

pm.table <- as.data.frame(modal.copy.number[(insecticide == 'PM') & (location %in% c('Baguida', 'Obuasi'))])
cat('\n\nAll populations PM:\n')
association.tests[['PM']][['combined']] <- glm.up(pm.table, genes.to.model, 'phenotype', control.for = 'location', glm.function = 'glmmTMB')
print(association.tests[['PM']][['combined']])

get.p.values <- function(results.object){
	mod <- results.object$final.model
	model.terms <- attributes(terms(mod))$term
	# If Coeae1f was in the final model, return it's marginal P-value from that model. Otherwise, 
	# return its P-value alone
	if ('Coeae1f' %in% model.terms)
		return(results.object$final.sig['Coeae1f', 'P'])
	else
		return(results.object$sig.alone['Coeae1f', 'P'])
}	

p.values.delta <- lapply(association.tests[['Delta']], get.p.values)
p.values.pm <- lapply(association.tests[['PM']], get.p.values)

# Calculate the odds ratios in Baguida, Obuasi and combined
get.odds.ratios <- function(results.object){
	mod <- results.object$final.model
	model.terms <- attributes(terms(mod))$term
	if (!('Coeae1f' %in% model.terms))
		return(c(NA, NA, NA))
	if (any(class(mod) == 'glmmTMB')){
		OR <- confint(mod)['Coeae1f', ] %>%
		      {1/exp(.)}
	}
	else {
		coeff <- mod$coefficients['Coeae1f']
		OR <- c(confint(mod)['Coeae1f', ],
		        coeff) %>%
		      {1/exp(.)}
	}
	# Have the mean in between the 2.5% and 97.5% intervals
	return(OR[c(1,3,2)])
}

odds.ratios.delta <- lapply(association.tests[['Delta']], get.odds.ratios)
odds.ratios.pm <- lapply(association.tests[['PM']], get.odds.ratios)

cat('\n\t######################')
cat('\n\t# Summary of results #')
cat('\n\t######################\n')

cat('\nSample sizes:\n')
print(population.popsize)

cat('\nCNV counts:\n')
population.modal.CNV.counts[, c('location', gene.name.conversion[Gene.name == 'Coeae1f', Gene.id])] %>%
print()

cat('\nCNV frequencies:\n')
print(CNV.frequencies)

cat('\nSignificance testing:\n')
cat('\n\tP-values:\n')
cat('\n\t\tDelta:\n')
print(p.values.delta)
cat('\n\t\tPM:\n')
print(p.values.pm)

cat('\n\tOdds ratios and confidence intervals (only reported for PM since delta was non-sig):\n')
cat('\n\t\tPM:\n')
print(odds.ratios.pm)

# Now show a contingency table of phenotype vs copy number
cat('\nPhenotype contingency table with Coeae1f in Baguida PM:\n')
Ba.PM.mcn <- modal.copy.number[list('Baguida', 'PM')]
Ba.coeae1f.table <- with(Ba.PM.mcn, table(Coeae1f,phenotype))
print(Ba.coeae1f.table)

cat('\nPhenotype contingency table with Coeae1f in Obuasi PM:\n')
Ob.PM.mcn <- modal.copy.number[list('Obuasi', 'PM')]
Ob.coeae1f.table <- with(Ob.PM.mcn, table(Coeae1f,phenotype))
print(Ob.coeae1f.table)

