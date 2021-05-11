#this file has 2 features:
#- repairing the java file new line characters (if a line starts with ., it will be apppended to previous line)
#- sorting the file alphabetically to make it easier for vimdiff to compare java files (only if mode=='sorted')
import sys

mode='sorted'

i=-1
previousline=''
results=[]
lines_list=open(sys.argv[1],'r').readlines()
lines=iter(lines_list)
for line in lines:
    print('current line',i,line)
    i+=1
    line_print=line.strip().strip('\n')
    if i<len(lines_list)-1:
        line_next=lines_list[i+1].strip().strip('\n')
        print('next line',line_next)
        if line_next.startswith('.'):
            print('next lien starts with .')
            line_print+=lines_list[i+1].strip().strip('\n')
            results.append(line_print)
            print('appending',line_print)
            i+=1
            next(lines)
        else:
            results.append(line_print)
    else:
        results.append(line_print)
    #thisline=line.strip()
    #if thisline.startswith('.'):
    #    results.append(previousline+thisline)
    #else:
    #    results.append(previousline)
    #previousline=line.strip()
if mode=='sorted':
    results=sorted(results)
print('Writing repaired java file to sorted.java')
with open('sorted.java','w') as of:
    for r in results:
        if len(r.strip())==0:
            continue
        of.write(r+'\n')
