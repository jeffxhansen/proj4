name = "Jeff"
name += " Hansen"
print(name[::-1])

print(len("ataagagtgattggcgtccgtacgtaccctttctactctcaaactcttgttagtttaaatctaatctaaactttataaacggcacttcctgtgtgtccat"))
MAXINDELS = 3
rows = 8
cols = 8
matrix = [[" " for j in range(cols)] for i in range(rows)]
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
        
print(matrix)
