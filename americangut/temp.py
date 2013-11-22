
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

def calculate_dimensions_bar(num_bars, bar_width = 0.5, axis_height = 3, 
    border = 0.1, title = 1, legend = 2, xlab = 0, ylab = 0, unit = 'inches'):
    """Determines the axis and figure dimensions for a bar chart.

    INPUTS:
        num_bars 

        border -- the size of white space to be added around the axis. 
                    DEFAULT: 0.1"

        title -- the height to add to the top of the figure for a title. This 
                    is separate from the border, which is added by default. 
                    DEFAULT: 1

        legend -- the width of the legend to the be added to the figure. 
                    DEFAULT: 2

        xlab -- the height to be added for labels along the x axis. DEFAULT: 0

        ylab -- the width to be added for labels along the y axis. DEFAULT: 0

        unit -- a string ('inches' or 'cm'), specifying the unit to be used 
                    in image generation. DEFAULT: inches

    OUTPUTS:
        axis_dimensions -- a Bbox class describing the axis position in the 
            figure

        figure_dimensions -- a 2 element tuple giving the width and height of 
            the figure in inches
    """

    # Preforms some sanity checks. 
    if num_bars < 1:
        raise ValueError ('There must be at least one group to plot.\nnum_bars'\
            ' must be a whole number greater than 1.')
    elif not round(num_bars) == num_bars:
        raise ValueError ('There cannot be partial categories. \nnum_bars must'\
            ' be a whole number greater than 1.')
 

   # Specifies a value for converting between units
    if unit == 'cm':
        conversion = float(1)/float(2.54)
    elif unit == 'inches':
        conversion = 1.0
    else:
        raise ValueError, 'unit must be "inches" or "cm".'

    # Determines the figure width
    axis_width = float(num_bars*bar_width)
    figure_width = float(border*2+ylab+axis_width+legend)
    figure_height = float(border*2+xlab+axis_height+title)

    figure_dimensions = (figure_width, figure_height)

    axis_left = float(ylab+border)/figure_width
    axis_right = float(ylab+border+axis_width)/figure_width
    axis_bottom = float(xlab+border)/figure_height
    axis_top = float(xlab+border+axis_height)/figure_height

    axis_dimensions = Bbox(array([[axis_left, axis_bottom],
                                  [axis_right, axis_top]]))

    return axis_dimensions, figure_dimensions

def render_bar_char(cat_table, cat_names, sample_names, axis_dims, 
    fig_dims, file_out = 'barchart', filetype = 'PDF', colormap = [[1, 1, 1]],
    legend = True, x_axis = True, x_tick_offset = 0.6, 
    bar_width = 0.8, x_min = -0.5, x_tick_interval = 1.0, y_lims = [0, 1], 
    y_tick_interval = 0.2, tick_font_size = 15, label_font_size = 20, 
    font_angle = 45, font_alignment = 'right'):
    """Creates a stacked barchart using the data in the category table.

    INPUTS:

    OUTPUT:
        The rendered figure is saved in the at the file_out location.       
    """

    # Preforms a sanity checks that the provided data is good
    (table_height, table_width) = cat_table.shape
    num_cats = len(cat_names)
    num_samples = len(sample_names)

    if not table_height == num_cats:
        raise ValueError ('The number of provided categories differ.')
    elif not table_width == num_samples:
        raise ValueError ('The number of samples differ.')

    # Sets up the x ticks 
    x_tick = arrange(0, num_samples)
    x_max = x_min + no_samples
    bar_left = x_tick - bar_width*0.5

    # Creates the figure
    fig = plt.figure(1, fig_dims)

    
    # Plots the data
    patches_watch = []
    for plot_count, category in enumerate(cat_table):
        bottom_bar = sum(category_table[:plot_count,:])
        faces = plt.bar(left = bar_left, 
                        height = category, 
                        width = bar_width, 
                        bottom = bottom_bar,
                        color = colormap[plot_count,:],
                        edgecolor = colormap[plot_count,:])
        patches_watch.append(faces[0])

    # Sets up the axis dimesnions and limits.
    ax1 = plt.gca()
    ax1.set_position(axis_dims)

     # The y-direction is reversed so the labels are in the same order as the 
    # colors in the legend
    plt.axis([x_min, x_max, y_lims[1], y_lims[0]])

     # Sets y axis labels
    y_tick_labels = (arange(y_max + y_tick_interval, y_min, \
        -y_tick_interval) - y_tick_interval)*100
    y_tick_labels[-1] = 0

    # Converts y label to text
    y_text_labels = [str(e) for e in y_tick_labels]

    y_tick_labels = ax1.set_yticklabels(y_text_labels, 
                                        size = tick_font_size)
    y_axis_label = ax1.set_ylabel('Frequency', 
                                  size = label_font_size)

    # Sets up the x-axis labels
    x_text_labels = sample_labels[:] # copy
    x_text_labels.insert(0, '') # insert at the head of the list
    
    x_tick_labels = ax1.set_xticklabels(x_text_labels, 
                                        size = tick_font_size,
                                        rotation = font_angle, 
                                        horizontalalignment = font_alignment)
    if legend:
        plt.figlegend(patches_watch, category_names, 'right')

    plt.savefig(file_out, format = file_format)


test_table = parse_biom_table(open('/Users/jwdebelius/Desktop/PythonTest1021/BODY_SITE_SPLIT/feces.biom'))
LEVEL = 2
[samp, cat, summ] = summarize_common_categories(test_table, LEVEL, limit = 1)

#print 'Sample Ids: %r' % samp
print sum(summ[:,:4])
# print 'Summary: %r' % summ




# def plot_stacked_phyla(taxonomy_table, taxonomy_headers, sample_labels, \
#   file_out, sample_ids=None, legend = True, x_axis=True):
#     """Creates a stacked taxonomy plot at the phylum level

#     INPUTS:
#         taxonomy_table -- a numpy array with sample information in the columns 
#                         and phylum frequency information in the rows

#         taxonomy_headers -- a  of the phyla which to the rows in the 
#                         taxonomy_table

#         sample_labels -- a list of the sample labels which correspond to the 
#                         rows in the taxonomy_table

#         file_out -- a string describing the filename for the output filename

#     OUTPUT:
#         The rendered figure is saved as a pdf in the at the file_out location.
#     """
#     # Colorbrewer python package. Can be loaded

#     # Colormap is taken from the colorbrewer    
#     COLORMAP = array([[0.8353, 0.2421, 0.3098],
#                       [0.9569, 0.4275, 0.2627],
#                       [0.9922, 0.6824, 0.3804],
#                       [0.9961, 0.8784, 0.5351],
#                       [0.9020, 0.9608, 0.5961],
#                       [0.6706, 0.8667, 0.6431],
#                       [0.4000, 0.7608, 0.6471],
#                       [0.1961, 0.5333, 0.7412],
#                       [0.3333, 0.3333, 0.3333]])

#     X_TICK_OFFSET = 0.6

#     BAR_WIDTH = 0.8

#     # Paramatizable the constants and the options for shape/size. User Passed size parameters

#     if legend and x_axis:
#         figure_dimensions = (8, 5)
#         axis_dimensions = Bbox(array([[0.2000, 0.3000],
#                                       [0.7000, 0.9000]]))
#     elif legend:
#         figure_dimensions = (8, 3.75)
#         axis_dimensions = Bbox(array([[0.2000, 0.1000],
#                                       [0.7000, 0.9000]]))
#     elif x_axis:
#         figure_dimensions = (6.2222, 5)
#         axis_dimensions = Bbox(array([[0.2000, 0.3000],
#                                       [0.9000, 0.9000]]))
#         #bah
#     else:
#         figure_dimensions = (6.22222, 3.75)
#         axis_dimensions = Bbox(array([[0.2000, 0.1000],
#                                       [0.9000, 0.9000]]))


#     X_MIN = -0.5
#     X_TICK_INTERVAL = 1.0

#     Y_MIN = 0
#     Y_MAX = 1.0
#     Y_TICK_INTERVAL = 0.2

#     TICK_FONT_SIZE = 15
#     LABEL_FONT_SIZE = 20

#     [no_phyla, no_samples] = taxonomy_table.shape

#     x_tick = arange(0,no_samples)
#     x_max = X_MIN+no_samples

#     sample_figure = plt.figure(1, figure_dimensions)


#     patches_watch = []

#     # Plots the data
#     for plot_count, phyla in enumerate(taxonomy_table):
#         already_added_index = arange(plot_count)
#         bottom_bar = sum(taxonomy_table[already_added_index,:])
#         faces = plt.bar(x_tick-BAR_WIDTH/2, phyla, BAR_WIDTH, \
#                         bottom = bottom_bar, color = COLORMAP[plot_count,:], \
#                         edgecolor=COLORMAP[plot_count,:])
#         patches_watch.append(faces[1])

#     # Sets up axis dimensions and limits
#     ax1 = plt.gca()
#     ax1.set_position(axis_dimensions)
    
#     # The y-direction is reversed so the labels are in the same order as the 
#     # colors in the legend
#     plt.axis([X_MIN, x_max, Y_MAX, Y_MIN])

#     # Sets y axis labels
#     y_tick_labels = (arange(Y_MAX + Y_TICK_INTERVAL, Y_MIN, \
#         -Y_TICK_INTERVAL) - Y_TICK_INTERVAL)*100
#     y_tick_labels[-1] = 0

#     # Converts y label to text
#     y_text_labels = [str(e) for e in y_tick_labels]

#     y_tick_labels = ax1.set_yticklabels(y_text_labels, size = TICK_FONT_SIZE)
#     y_axis_label = ax1.set_ylabel('Frequency', size = LABEL_FONT_SIZE)

#     # Sets up the x-axis labels
#     x_text_labels = sample_labels[:] # copy
#     x_text_labels.insert(0, '') # insert at the head of the list
    
#     x_tick_labels = ax1.set_xticklabels(x_text_labels, size = TICK_FONT_SIZE,\
#         rotation = 45, horizontalalignment = 'right')

#     # Adds the legend
#     if legend:
#         plt.figlegend(patches_watch, taxonomy_headers, 'right')


#     plt.savefig(file_out, format = 'pdf')


