in_f = open('all_dem.txt')
orig_list = open('left_overs.txt')
ou_f = open('all.cnv')
map_match = {}
map_test = [] 
for line in orig_list:
    map_test.append(line.strip())
for i, line in enumerate(in_f):
    map_match[(line.strip())] = map_test[i]

find_li=[]
for line in ou_f:
    find_li.append(line.strip()) 

for line in map_match.keys():
    line=line.strip()
    if any(line in s for s in find_li):
        print(map_match[line])
