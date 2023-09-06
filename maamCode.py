import numpy as np
import random
import matplotlib.pyplot as plt
import cv2
print("Generating random value for V between 0 and 1 ")
del_t = 0.01
L=20.0
# L is the dimension of the graph
v = 1.0
print("V = "+str(v)+" generated")
k = 1.0
arr = []
n_part = 100
# npart is the number of molecules present on the surface
n_iter = int(input("Enter the number of iterations :"))
vecx = [[] for i in range(n_part)]
vecy = [[] for i in range(n_part)]
for i in range(n_part):
    x = random.random()*L
    y = random.random()*L
    vecx[i].append(x)
    vecy[i].append(y)
    arr.append((x,y))
# above function generates the arr of molecules in the graph
endo_rate = 10
exo_rate = 10
sigma_rate = 5
freq = 100
recycleFreq = 30
recycled_particles = 0
recycleRate = 0.3
l0=1.1 # touching distance of molecules
img_names = []
for t in range(n_iter):
    # print(len(arr))
    # print("vecx is "+str(vecx))
    # print("vecy is "+str(vecy))
    if t%freq==0:
        # this code runs every 100th iterations!
        # basically saves the current configuration of the surface
    # print("vecx is "+str(vecx))
    # print("vecy is "+str(vecy))
        fig, ax = plt.subplots()
        to_printx =[]
        to_printy =[]
        for i in range(n_part):
            to_printx.append(vecx[i][-1])
            to_printy.append(vecy[i][-1])
        ax.plot(to_printx, to_printy,'o')
        fig.savefig(f"graph_{t}.png")
        plt.close(fig)
        img_names.append(f"graph_{t}.png")

 # removing molecules
        nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
        print(nd)
        while(nd>=n_part):
            nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
        recycled_particles = int(np.round((n_part - len(arr))*recycleRate))
        nd -= recycled_particles
        # this above code used to essentially make the nd (number of removed molecules less than n_part)
        removed_ind = set()
        while(len(removed_ind)<nd):
            x = random.randint(0, n_part-1)
            removed_ind.add(x)
        # removed_ind has now all the indexes that are to be removed, but error here as x can be repeated
        new_indices = []
        for i in range(n_part):
            if i not in removed_ind:
                new_indices.append(i)
        # above comment got taken care of here, but still if there was a repeat the nd would be more than the actual removed molecules
        tvecx=[]
        tvecy = []
        tarr=[]
        for ind in new_indices:
            tvecx.append(vecx[ind])
            tvecy.append(vecy[ind])
            tarr.append(arr[ind])
        vecx = tvecx
        vecy = tvecy
        arr=tarr
        print("after removing molecules")
        print(len(arr))
        
        
        # replacing the actual list and positions after removing molecules
        # adding molecules
        na = int(np.round(random.gauss(exo_rate,sigma_rate))) if exo_rate != 0 else 0
        for i in range(na):
            x = random.random()*L
            y = random.random()*L
            vecx.append([x])
            vecy.append([y])
            arr.append((x,y))
        n_part = len(vecx)
        print("after adding molecules")
        print(len(arr))
        # molecules are now added, so n_part is increased

    # adding recyled particles
    elif (t-recycleFreq)%freq==0:
        # assuming constant endocytosis rate
        for i in range(recycled_particles):
            x = random.random()*L
            y = random.random()*L
            vecx.append([x])
            vecy.append([y])
            arr.append((x,y))
        n_part = len(vecx)
        # exo_rate = exo_rate - recycled_particles if exo_rate - recycled_particles >= 0 else 0
    fxij = np.zeros(n_part)
    fyij = np.zeros(n_part)
    # print(n_part)
    # print("size of vecx is "+str(len(vecx)))
    # print("size of arr is "+str(len(arr)))
    # print()
    for i in range(n_part-1):
        for j in range(i+1,n_part):
            (xi,yi)=arr[i]
            (xj,yj)=arr[j]
            delxij = (xi-xj)
            delyij = (yi-yj)
            delxij = delxij - np.round(delxij/L)*L
            delyij = delyij - np.round(delyij/L)*L
            lij = np.sqrt(delxij**2.0 + delyij**2.0)
            if lij <= 2.0*l0:
                fij = -k*(lij-l0)
                fxij[i] += fij*delxij/abs(lij)
                fxij[j] -= fij*delxij/abs(lij)
                fyij[i] += fij*delyij/abs(lij)
                fyij[j] -= fij*delyij/abs(lij)
    
    for i in range(n_part):
        (x,y)=arr[i]
        valid = False
        while(not valid):
            temp = random.random()
            theta = 2*np.pi*temp
            x = x + (fxij[i] + v * np.cos(theta)) * del_t
            y = y + (fyij[i] + v * np.sin(theta)) * del_t
            # if(x<0):
            # x+=L
            # if(x>L):
            # x-=L
            # if(y<0):
            # y+=L
            # if(y>L):
            # y-=L
            temp=True
            valid=True
        vecx[i].append(x)
        vecy[i].append(y)
        arr[i] = (x, y)

print("pre-processing done!")
# Set up video codec and output file name
codec = cv2.VideoWriter_fourcc(*"mp4v")
out = cv2.VideoWriter("output.mp4", codec, 18, (640, 480))
for i in range(len(img_names)):
    # Read the saved image file and add it to the video
    img = cv2.imread(img_names[i])
    out.write(img)
plt.show()
# Release the video writer and cleanup
out.release()
cv2.destroyAllWindows()
print("done!")