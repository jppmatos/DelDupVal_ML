#python seg_dup.py segmental_duplications_file > outfile
from sys import argv


def read_dups(infile):
        f=open(infile)
        l=[]
        dic={}
        a=1
        for i in f:
                line=i.split("\t")
                if [line[1], line[2], line[3]] not in l and [line[7], line[8], line[9]] not in l:
                        dic["seg_dup"+str(a)]=[[line[1], line[2], line[3]], [line[7], line[8], line[9]]]
                        a+=1
                        l.append([line[1], line[2], line[3]])
                        l.append([line[7], line[8], line[9]])
        f.close()
        return dic

def write_dic(dic):
        #c=0
        for key,value in dic.items():
                print(key+"\t"+"\t".join(value[0])+"\t"+"\t".join(value[1]))
                c+=1
                #if  c > 6:
                    #break
write_dic(read_dups(argv[1]))

