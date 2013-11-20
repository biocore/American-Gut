
from biom.parse import parse_biom_table
from numpy import mean, array

def summarize_common_categories(otu_table, level, metadata_category = 'taxonomy',
	limit_mode = 'COMPOSITE', limit = 1.0):
	
	"""Identifies the most common taxa in a population using variable limits

	This method uses a composite score to account for both the average frequency
	and the percentage of the population in which the sample is present. This 
	hopes to be more robust than either the average or fraction present alone.

	INPUTS:
		biom_table -- a sparse biom table to be summarized

		metadata_category -- a description of the metadata category over which 
					the data will be summarized.

        level -- an integer corresponding to the taxonomic level (or other meta 
                    data category level) at which data should be summarized.
                    The argument will only take a single integer.

        limit_mode -- a string describing how the scoring limit which should be 
        			used. Options are 'COMPOSITE', 'AVERAGE', 'COUNTS' or 
        			'NUMBER', 'NONE'. 
        			COMPOSITE keeps samples whose composite scores greater than 
        				the limit. (Default limit is 1)
        			AVERAGE keeps samples whose average are greater than the 
        				specified limit. (Default limit is 0.01)
        			COUNTS keep groups which are represented in at least the 
        				specified fraction of sample (default limit is 0.01)
        			NONE sorts groups alphabetically and keeps all.        			

        limit -- the numeric lower limit for common taxa

        logging -- a binary value specifying if a table of the frequency, 
        			composite value, and metadata category should be printed and
        			retained.

        log_file -- the filename where the log file should be saved.

    OUTPUTS
        sample_ids -- a list of the sample ids contained in the matrix 
        			(column header)
    	
    	common_categories -- a list of the common categories in the (row header)

    	category_summary -- a numpy array of the category values.

    	scores -- a sorted list of all the taxa and their scores

	"""	
	# Sets up positions
	COMPOSITE_CONSTANT = 100^2
	header = ['%s', 'Composite Score', 'Average Freq', 'Population Freq']

	if limit_mode == 'COMPOSITE':
		score_position = 3
		method = 'using the composite score'
		limit_string = 'with %g as a minimum value.' % limit
	elif limit_mode == 'AVERAGE':
		score_position == 1
		method = 'using the average value'
		limit_string = 'with %g as a minimum value.' % limit
	elif limit_mode == 'COUNTS':
		score_position = 2
		method = 'using the population frequency'
		limit_string = 'with %g as a minimum value.' % limit
	elif limit_mode == 'NONE':
		score_position = 0
		method = 'by sorting alphabetically'
		limit_string = 'without a minimum value.'
	else:
		raise ValueError, 'limit_mode is not a supported option.'

	# Prealocates
	scoring_all = []
	common_categories = []
	category_summary = []
	table_row = []
	other = []

	# Gets the Sample IDs
	sample_ids = list(otu_table.SampleIds)
	num_samples = len(sample_ids)


	# Normalizes the data by the number of observations so relative frequencies 
	# are used.
	otu_table_norm = otu_table.normObservationBySample()
	otu_table_norm.__class__
	# Collapses the OTUs into category summaries using the correct levels
	for (bin, table) in otu_table_norm.binObservationsByMetadata(lambda x: 
			x[metadata_category][:level]):		
		# Pulls out the sample data for the group
		group_value = table.sum('sample')
		group_binary = group_value > 0

		# Calculates presence scores
		average_freq = mean(group_value)
		fraction_pres = float(sum(group_binary))/float(num_samples)
		composite = average_freq*fraction_pres*COMPOSITE_CONSTANT
		score_row = [bin, average_freq, fraction_pres, composite, group_value]		

		# Adds the scores to the watch matrix
		if fraction_pres > 0:
			scoring_all.append(score_row)

	# Sorts based on scoring method
	scores = sorted(scoring_all, key = lambda score: scoring_all[score_position])

	# Identifies rows which meet the scoring criteria
	for score in scores:
		#print score[0:3]
		if score[score_position] > limit:
			common_categories.append(score[0])
			category_summary.append(score[4])

	common_categories.append((u'Other'))
	category_summary.append(1-sum(category_summary))

	# Logs the data if necessary
	if len(common_categories) == 0:
		raise ValueError, ('Limit too high! No common categories could be'\
			' identified.')


	# Returns the values
	return sample_ids, common_categories, array(category_summary)




test_table = parse_biom_table(open('/Users/jwdebelius/Desktop/PythonTest1021/BODY_SITE_SPLIT/feces.biom'))
LEVEL = 2
[samp, cat, summ] = summarize_common_categories(test_table, LEVEL, limit = 1)

#print 'Sample Ids: %r' % samp
print sum(summ[:,:4])
# print 'Summary: %r' % summ

