import numpy as np
import matplotlib.pyplot as plt


def file_open(file):
    import sys
    data = []
    try:
        f = open(file, 'r', encoding='utf-8')
    except Exception:
        print("open error. not found file:", str(file))
        sys.exit(1)
    for line in f:
        line = line.strip() #前後空白削除
        line = line.replace('\n','') #末尾の\nの削除
        line = line.split(",") #分割
        data.append(line)
    f.close()
    return data

#左から順番に、v,a,E,ri,ru,lu,med
data=file_open('log/logfile5.txt')

a=[1.8]
v1,v2,v3,v4,v5,v6,v7,v8,v9=[],[],[],[],[],[],[],[],[]
v=[v1,v2,v3,v4,v5,v6,v7,v8,v9]
for i in range(len(data)):
    for j in range(len(a)):
        if data[i][1]==str(a[j]):
            v[j].append([round(float(data[i][0]),4),
                         round(float(data[i][2]),4),
                         round(float(data[i][3]),4),
                         round(float(data[i][4]),4),
                         round(float(data[i][5]),4),
                         round(float(data[i][6]),4),
                         round(float(data[i][7]),4)])

for j in range(len(v)):
    v[j].sort(reverse=False, key=lambda x:x[0])
#uは左から順番に、v,E,ri,ru,lu,med

# f=open("LOGFILE.txt","w")
# for i in range(len(v)):
#     if len(v)==0:break
#     for j in range(len(v[i])):
#         f.write("{},{},{}\n".format( v[i][j][0],a[i],v[i][j][5] ))
# f.close()


# for i in range(len(v)):
#     u=np.array(v[i])
#     if len(u)!=0:

#         fig=plt.figure(figsize=(16,6))

#         ax1=fig.add_subplot(2,3,1)
#         ax1.scatter(u[:,0],u[:,1],s=15)
#         ax1.set_title("{}".format(a[i]))

#         ax2=fig.add_subplot(2,3,2)
#         ax2.scatter(u[:,0],u[:,5],s=15)

#         ax3=fig.add_subplot(2,3,3)
#         ax3.scatter(u[:,0],u[:,2],s=15)

#         ax4=fig.add_subplot(2,3,4)
#         ax3.scatter(u[:,0],u[:,3],s=15)

#         ax5=fig.add_subplot(2,3,5)
#         ax3.scatter(u[:,0],u[:,4],s=15)

#         ax6=fig.add_subplot(2,3,6)
#         ax6.scatter(u[:,0],u[:,6],s=15)
        
#         plt.show()



u=np.array(v[0])
plt.figure(figsize=(10,6))
plt.rcParams["font.size"] = 20

# plt.scatter(u[:,0],u[:,1],color='blue',label="$k_1$",s=35)
plt.scatter(u[:,0],u[:,2],label="$k_1$",s=35)
plt.scatter(u[:,0],u[:,3],label="$k_2$",s=35)
plt.scatter(u[:,0],u[:,4],label="$k_3$",s=35)
plt.legend()

plt.xlabel("$V_0$",fontsize=20)
# plt.xticks([0,0.3,0.6,0.9,1.2,1.5])
# plt.yticks([0,0.2,0.4,0.6,0.8])
plt.xticks([0,0.5,1.0,1.5,2.0,2.5,3,3.5,4])
plt.yticks([0,0.2,0.4,0.6,0.8,1.0])
plt.ylabel("$|S(k)|$",fontsize=20)
# plt.ylabel("$E_{period}/E_{\Omega=0}$",fontsize=20)

min=0
max=1
plt.ylim(min,max)
plt.vlines(x=2.35, ymin=min, ymax=max,colors="black",linestyle='dashed', linewidth=1)
plt.tight_layout()
plt.show()
    
