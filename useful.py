import time

# Used to compute the bandwidth for banded version
MAXINDELS = 3

# Used to implement Needleman-Wunsch scoring
MATCH = -3
INDEL = 5
SUB = 1

INS_L = 'L'
DEL_U = 'U'
SUB_D = 'D'


class Edit_Distancizer:

	def __init__(self):
		self.edit_distance = []
		self.prev_matrix = []

	def initialize_matrices(self, seq1, seq2, banded):
		# get rows and columns
		rows = len(seq1)+1
		cols = len(seq2)+1

		if self.edit_distance != [] and self.prev_matrix != []:

			width = len(self.edit_distance[0])
			height = len(self.edit_distance)

			if width == cols and height == rows:
				return

		# create the computation matrix as a 2D array
		# create the previous matrix as a 2d array
		if banded:
			fill = float('inf')
		else:
			fill = 0
		self.edit_distance = [[fill
								for j in range(cols)] for i in range(rows)]
		self.prev_matrix = [['' for j in range(cols)] for i in range(rows)]

	def best_option(self, i_char, j_char, i, j):
		ins_option = self.edit_distance[i][j-1] + INDEL
		del_option = self.edit_distance[i-1][j] + INDEL
		if i_char == j_char:
			sub_option = self.edit_distance[i-1][j-1] + MATCH
		else:
			sub_option = self.edit_distance[i-1][j-1] + SUB

		dist = ins_option
		prev = INS_L
		if ins_option < dist:
			dist = ins_option
			prev = INS_L
		if del_option < dist:
			dist = del_option
			prev = DEL_U
		if sub_option < dist:
			dist = sub_option
			prev = SUB_D

		return dist, prev

	def get_alignments(self, seq1, seq2):

		i = len(self.prev_matrix)-1
		j = len(self.prev_matrix[0])-1
		pointer = self.prev_matrix[i][j]
		align1 = ""
		align2 = ""
		while pointer != "":
			if pointer == SUB_D:
				align1 += seq1[i-1]
				align2 += seq2[j-1]
				i -= 1
				j -= 1
			if pointer == INS_L:
				align1 += "-"
				align2 += seq2[j-1]
				j -= 1
			if pointer == DEL_U:
				align1 += seq1[i-1]
				align2 += "-"
				i -= 1
			pointer = self.prev_matrix[i][j]

		# reverse the alignments so they point front to back
		align1 = align1[::-1]
		align2 = align2[::-1]

		return align1, align2

	def alignments_score_banded(self, seq1, seq2):
		# get rows and columns
		rows = len(seq1)+1
		cols = len(seq2)+1

		# create the computation matrix as a 2D array
		# create the previous matrix as a 2d array
		t1 = time.time()
		self.initialize_matrices(seq1, seq2, banded=True)
		t2 = time.time()
		print("{} seconds to initialize matrices".format(t2-t1))

		# initialize the top row and left column sequentially
		# inititialize the previous matrix in similar fashion
		for i in range(MAXINDELS+1):
			self.edit_distance[i][0] = i*INDEL
			self.prev_matrix[i][0] = DEL_U
		for j in range(MAXINDELS+1):
			self.edit_distance[0][j] = j*INDEL
			self.prev_matrix[0][j] = INS_L
		self.prev_matrix[0][0] = ''

		# iterate through the matrix updating each cell in
		# computation matrix and prev matrix
		t3 = time.time()
		for i in range(1, rows):
			start = i - MAXINDELS
			end = i + MAXINDELS + 1
			if start < 1:
				start = 1
			if end > rows:
				end = rows
			for j in range(start, end):
				i_char = seq1[i-1]
				j_char = seq2[j-1]
				dist, prev = self.best_option(i_char, j_char, i, j)
				self.edit_distance[i][j] = dist
				self.prev_matrix[i][j] = prev
		t4 = time.time()
		print("{} seconds to run the banded algorithm".format(t4-t3))
		print()

		print(self.edit_distance[-1][-1])

		# get the alignments
		align1, align2 = self.get_alignments(seq1, seq2)

		# return alignment1, alignment2, and
		# score (bottom right value in computation table)
		return align1[:100], align2[:100], self.edit_distance[-1][-1]

	def alignments_score_unbanded(self, seq1, seq2):

		# get rows and columns
		rows = len(seq1)+1
		cols = len(seq2)+1

		# create the computation matrix as a 2D array
		# create the previous matrix as a 2d array
		t1 = time.time()
		self.initialize_matrices(seq1, seq2, banded=False)
		t2 = time.time()
		print("{} seconds to initialize matrices".format(t2-t1))

		# initialize the top row and left column sequentially
		# inititialize the previous matrix in similar fashion
		for i in range(rows):
			self.edit_distance[i][0] = i*INDEL
			self.prev_matrix[i][0] = DEL_U
		for j in range(cols):
			self.edit_distance[0][j] = j*INDEL
			self.prev_matrix[0][j] = INS_L
		self.prev_matrix[0][0] = ''

		# iterate through the matrix updating each cell in
		# computation matrix and prev matrix
		t3 = time.time()
		for i in range(1, rows):
			for j in range(1, cols):
				i_char = seq1[i-1]
				j_char = seq2[j-1]
				dist, prev = self.best_option(i_char, j_char, i, j)
				self.edit_distance[i][j] = dist
				self.prev_matrix[i][j] = prev
		t4 = time.time()
		print("{} seconds to run unbanded".format(t4-t3))

		print(self.edit_distance[-1][-1])

		# get the alignments
		align1, align2 = self.get_alignments(seq1, seq2)

		# return alignment1, alignment2, and
		# score (bottom right value in computation table)
		return align1[:100], align2[:100], self.edit_distance[-1][-1]


# This is the method called by the GUI.  _seq1_ and _seq2_ are two sequences to be aligned, _banded_ is a boolean that tells
# you whether you should compute a banded alignment or full alignment, and _align_length_ tells you
# how many base pairs to use in computing the alignment


	def align(self, seq1, seq2, banded):
		self.banded = banded

		if not banded:
			alignment1, alignment2, score = \
				self.alignments_score_unbanded(seq1, seq2)
		else:
			alignment1, alignment2, score = \
				self.alignments_score_banded(seq1, seq2)
    
		print(alignment1)
		print(alignment2)
		for line in self.edit_distance:
			print(line)
		for line in self.prev_matrix:
			print(line)

		return {'align_cost': score, 'seq1': alignment1, 'seq2': alignment2, "dist_matrix":self.edit_distance, "prev_matrix":self.prev_matrix}


ed = Edit_Distancizer()

results = ed.align("ATGCC", "TACGCA", False)