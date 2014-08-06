SILENT_MODE
BLOCK_FILE jp_3981BYL.concise.blc
MAX_NSEQ 362
MAX_INPUT_LEN 364
OUTPUT_FILE jp_3981BYL.concise.ps
PORTRAIT
POINTSIZE 8
IDENT_WIDTH 12
X_OFFSET 2
Y_OFFSET 2
DEFINE_FONT 0 Helvetica		DEFAULT 
DEFINE_FONT 1 Helvetica		REL		0.75   
DEFINE_FONT 7 Helvetica		REL		0.6
DEFINE_FONT 3 Helvetica-Bold	DEFAULT    
DEFINE_FONT 4 Times-Bold    	DEFAULT   
DEFINE_FONT 5 Helvetica-BoldOblique	DEFAULT 
#
DEFINE_COLOUR 3  1 0.62 0.67	# Turquiose
DEFINE_COLOUR 4  1 1 0		# Yellow
DEFINE_COLOUR 5  1 0 0		# Red
DEFINE_COLOUR 7  1 0 1		# Purple
DEFINE_COLOUR 8  0 0 1		# Blue
DEFINE_COLOUR 9  0 1 0		# Green
DEFINE_COLOUR 10 0.41 0.64 1.00	# Pale blue 
DEFINE_COLOUR 11 0.41 0.82 0.67	# Pale green 
DEFINE_COLOUR 50 0.69 0.18 0.37	# Pink (helix)
DEFINE_COLOUR 51 1.00 0.89 0.00	# Gold (strand)
NUMBER_INT 10
SETUP
#
# Highlight specific residues.
# Avoid highlighting Lupas 'C' predictions by
# limiting the highlighting to the alignments 
Scol_CHARS	C 1 1 175 351   4
Ccol_CHARS	H ALL    5
Ccol_CHARS	P ALL    8
SURROUND_CHARS LIV   ALL
#
# Replace known structure types with whitespace
SUB_CHARS 1 352 175 361 H SPACE
SUB_CHARS 1 352 175 361 E SPACE
SUB_CHARS 1 352 175 361 - SPACE
STRAND 53 355 55
COLOUR_TEXT_REGION 53 355 55 355 51
STRAND 92 355 95
COLOUR_TEXT_REGION 92 355 95 355 51
HELIX 58 355 70
COLOUR_TEXT_REGION 58 355 70 355 50
HELIX 76 355 86
COLOUR_TEXT_REGION 76 355 86 355 50
HELIX 101 355 108
COLOUR_TEXT_REGION 101 355 108 355 50
STRAND 52 360 55
COLOUR_TEXT_REGION 52 360 55 360 51
STRAND 91 360 97
COLOUR_TEXT_REGION 91 360 97 360 51
HELIX 58 360 70
COLOUR_TEXT_REGION 58 360 70 360 50
HELIX 76 360 86
COLOUR_TEXT_REGION 76 360 86 360 50
HELIX 100 360 107
COLOUR_TEXT_REGION 100 360 107 360 50
STRAND 54 361 54
COLOUR_TEXT_REGION 54 361 54 361 51
HELIX 58 361 70
COLOUR_TEXT_REGION 58 361 70 361 50
HELIX 77 361 85
COLOUR_TEXT_REGION 77 361 85 361 50
HELIX 92 361 94
COLOUR_TEXT_REGION 92 361 94 361 50
HELIX 103 361 109
COLOUR_TEXT_REGION 103 361 109 361 50