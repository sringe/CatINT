
def replace_line(file_name, line_num, text):
    lines = open(file_name, 'r').readlines()
    lines[line_num] = text
    out = open(file_name, 'w')
    out.writelines(lines)
    out.close()

def insert_line(file_name,line_num,text):
    f = open(file_name, "r")
    contents = f.readlines()
    f.close()
    contents.insert(line_num, text)
    f = open(file_name, "w")
    contents = "".join(contents)
    f.write(contents)
    f.close() 
