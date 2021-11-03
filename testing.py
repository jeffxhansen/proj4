'''
MAXINDELS = 3
seq1 = "ACAATCC"
seq2 = "AGCATGC"
rows = len(seq1) + 1
cols = 2*MAXINDELS+1
matrix = [["-" for j in range(cols)] for i in range(rows)]

for i in range(MAXINDELS+1):
    matrix[i][0] = "_"
for j in range(MAXINDELS+1):
    if j ==0:
        matrix[0][j] = "_"
    else:
        matrix[0][j] = seq2[j-1]

k = 2*MAXINDELS+1
offset = -(MAXINDELS)
for i in range(1, rows):
    if offset < 0:
        start = 1
        end = k + offset
    elif offset + k:
        start = 0
        end = k
    else:
        start = 0
        end = rows - i + MAXINDELS
    print(start, end, offset)
    for j in range(start, end+1):
        if offset < 0:
            index = j-1
        else:
            index = j + offset
        if index == 0:
            ch = "_"
        else:
            ch = seq2[index-1]
        matrix[i][index] = ch
        
    offset += 1

for line in matrix:
    print("".join(line))
'''

MAXINDELS = 3
rows = 10
cols = 10
matrix = [["-" for j in range(cols)] for i in range(rows)]
for i in range(MAXINDELS+1):
    matrix[i][0] = "*"
for j in range(MAXINDELS+1):
    matrix[0][j] = "*"

for i in range(1, rows):
    start = i - 3
    end = i + 4
    if start < 1:
        start = 1
    if end > rows:
        end = rows
    for j in range(start, end):
        matrix[i][j] = "*"
        
for line in matrix:
    print("".join(line))
    
print()
    
MAXINDELS = 3
k = 2*MAXINDELS+1
seq1 = "ACAATCC"*2
seq2 = "AGCATGC"*2
rows = len(seq1)+1
cols = k
matrix = [["-" for j in range(cols)] for i in range(rows)]
for i in range(MAXINDELS+1):
    matrix[i][0] = "*"
for j in range(MAXINDELS+1):
    matrix[0][j] = "*"

offset = -(MAXINDELS)
offsets = ["(NA, NA, NA)"]

for i in range(1, rows):
    start = i - MAXINDELS
    end = i + MAXINDELS + 1
    if start < 1:
        start = 1
    if end > rows:
        end = rows
    '''
    if offset > 0:
        start -= offset
        end -= offset
    '''
    if offset > -1:
        start -= (offset+1)
        end -= (offset+1)
        frontChange = k-end
        
        if end < k:
            start = frontChange
            end = k
        
    
    for j in range(start, end):
        matrix[i][j] = "*"
    offsets.append((offset, start, end))
    offset += 1

for i, line in enumerate(matrix):
    print(str(offsets[i]) + ":\t" + "".join(line))

print()
