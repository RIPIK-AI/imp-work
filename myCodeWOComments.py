import numpy as np
import random
import matplotlib.pyplot as plt
import cv2
from datetime import datetime
def Rab4Recycle(removed_particles_entering_rab4, rab4_efficiency):
    return int(np.round(removed_particles_entering_rab4*rab4_efficiency))
def TubularRecycle(removed_particles_entering_tubular, tubular_efficiency):
    return int(np.round(removed_particles_entering_tubular*tubular_efficiency))
def Rab4Efficiency(n_part):
    if n_part < 40:
        return 0.5
    elif n_part > 100:
        return 0.2
    else:
        return -0.005 * n_part + 0.7
def wildType(n_iter, mechanism_type):
    currentDate = datetime.now().date()
    currentTime = datetime.now().time()
    currentTime = str(currentTime)[:8]
    part_arr = []
    numbers = []
    L=20.0
    v = 1.0
    k = 1.0
    arr = []
    n_part = 100
    vecx = [[] for i in range(n_part)]
    vecy = [[] for i in range(n_part)]
    for i in range(n_part):
        x = random.random()*L
        y = random.random()*L
        vecx[i].append(x)
        vecy[i].append(y)
        arr.append((x,y))
    endo_rate = 10
    exo_rate = 10
    sigma_rate = 5
    freq = 100
    rab4_efficiency = 0.2
    lamRecycleFreqRab4 = 20
    lamRecycleFreqTubular = 40
    recycleFreqTubular = 40 
    recycleFreqRab4 = 20
    tubular_efficiency = 0.4
    img_names = []
    for t in range(n_iter):
        if t%freq==0:
            fig, ax = plt.subplots()
            to_printx =[]
            to_printy =[]
            for i in range(n_part):
                to_printx.append(vecx[i][-1])
                to_printy.append(vecy[i][-1])
            ax.plot(to_printx, to_printy,'o')
            fig.savefig(f"images_my/{mechanism_type}/graph_{t}.png")
            plt.close(fig)
            img_names.append(f"images_my/{mechanism_type}/graph_{t}.png")
            nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
            while(nd>=n_part):
                nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
            recycleFreqTubular = np.random.poisson(lam=lamRecycleFreqTubular)
            recycleFreqRab4 = np.random.poisson(lam=lamRecycleFreqRab4)
            recycled_particles_rab4 = Rab4Recycle(removed_particles_entering_rab4= nd*0.5, rab4_efficiency=Rab4Efficiency(n_part))
            recycled_particles_tubular = TubularRecycle(removed_particles_entering_tubular=nd*0.5, tubular_efficiency=tubular_efficiency)
            removed_ind = set()
            while(len(removed_ind)<nd):
                x = random.randint(0, n_part-1)
                removed_ind.add(x)
            new_indices = []
            for i in range(n_part):
                if i not in removed_ind:
                    new_indices.append(i)
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
            na = int(np.round(random.gauss(exo_rate,sigma_rate))) - (recycled_particles_rab4 + recycled_particles_tubular)
            for i in range(na):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
            part_arr.append(n_part)
        elif (t+recycleFreqRab4)%freq==0:
            for i in range(recycled_particles_rab4):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
        elif (t+recycleFreqTubular)%freq==0:
            for i in range(recycled_particles_tubular):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
        if t%10 == 0:
            numbers.append(n_part)
    total_sum = 0
    for part in range(0, len(part_arr)):
        total_sum += part_arr[part]
    print(total_sum/len(part_arr))
    plt.figure(figsize=(10,10))
    plt.plot(numbers)
    plt.ylim(bottom=0)
    plt.title('Line Graph of E-cad molecules on surface vs Time')
    plt.xlabel('Time (Every 10th iteration)')
    plt.ylabel('Number of Particles')
    plt.savefig(f"graphs/{mechanism_type}/graph_{currentDate} {currentTime}.png")
    return img_names
def Rab11(n_iter, mechanism_type):
    currentDate = datetime.now().date()
    currentTime = datetime.now().time()
    currentTime = str(currentTime)[:8]
    part_arr = []
    numbers = []
    L=20.0
    arr = []
    n_part = 100
    vecx = [[] for i in range(n_part)]
    vecy = [[] for i in range(n_part)]
    for i in range(n_part):
        x = random.random()*L
        y = random.random()*L
        vecx[i].append(x)
        vecy[i].append(y)
        arr.append((x,y))
    endo_rate = 10
    exo_rate = 10
    sigma_rate = 5
    freq = 100
    recycleFreq = 30
    lamRecycleFreqRab4 = 20
    lamRecycleFreqTubular = 40
    recycleFreqTubular = 40 
    recycleFreqRab4 = 20
    rab4_efficiency = 0.2
    tubular_efficiency = 0.4
    img_names = []  
    for t in range(n_iter):
        if t%freq==0:
            fig, ax = plt.subplots()
            to_printx =[]
            to_printy =[]
            for i in range(n_part):
                to_printx.append(vecx[i][-1])
                to_printy.append(vecy[i][-1])
            ax.plot(to_printx, to_printy,'o')
            fig.savefig(f"images_my/{mechanism_type}/graph_{t}.png")
            plt.close(fig)
            img_names.append(f"images_my/{mechanism_type}/graph_{t}.png")
            nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
            while(nd>=n_part):
                nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
            recycleFreqTubular = np.random.poisson(lam=lamRecycleFreqTubular)
            recycleFreqRab4 = np.random.poisson(lam=lamRecycleFreqRab4)
            recycled_particles_rab4 = Rab4Recycle(removed_particles_entering_rab4= nd*0.5, rab4_efficiency=Rab4Efficiency(n_part))
            recycled_particles_rab4_supposed = Rab4Recycle(removed_particles_entering_rab4= nd*0.5, rab4_efficiency=0.2)
            recycled_particles_tubular = TubularRecycle(removed_particles_entering_tubular=nd*0.5, tubular_efficiency=tubular_efficiency)
            recycled_particles_tubular_supposed = TubularRecycle(removed_particles_entering_tubular=nd*0.5, tubular_efficiency=0.4)
            removed_ind = set()
            while(len(removed_ind)<nd):
                x = random.randint(0, n_part-1)
                removed_ind.add(x)
            new_indices = []
            for i in range(n_part):
                if i not in removed_ind:
                    new_indices.append(i)
            tvecx=[]
            tvecy = []
            tarr=[]
            if len(new_indices) <= 50:
                rab4_efficiency = 0.3
            elif len(new_indices) >= 75:
                rab4_efficiency = 0.1
            for ind in new_indices:
                tvecx.append(vecx[ind])
                tvecy.append(vecy[ind])
                tarr.append(arr[ind])
            vecx = tvecx
            vecy = tvecy
            arr=tarr
            na = int(np.round(random.gauss(exo_rate,sigma_rate))) - (recycled_particles_rab4_supposed + recycled_particles_tubular_supposed)
            for i in range(na):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
            part_arr.append(n_part)
        elif (t+recycleFreqRab4)%freq==0:
            for i in range(recycled_particles_rab4):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
        if t%100 == 0:
            numbers.append(n_part)
    total_sum = 0
    for part in range(0, len(part_arr)):
        total_sum += part_arr[part]
    print(total_sum/len(part_arr))
    plt.figure(figsize=(10,10))
    plt.plot(numbers)
    plt.ylim(bottom=0)
    plt.title('Line Graph of E-cad molecules on surface vs Time')
    plt.xlabel('Time (Every 10th iteration)')
    plt.ylabel('Number of Particles')
    plt.savefig(f"graphs/{mechanism_type}/graph_{currentDate} {currentTime}.png")
    return img_names
def RabX1(n_iter, mechanism_type):
    currentDate = datetime.now().date()
    currentTime = datetime.now().time()
    currentTime = str(currentTime)[:8]
    part_arr = []
    numbers = []
    L=20.0
    arr = []
    n_part = 100
    vecx = [[] for i in range(n_part)]
    vecy = [[] for i in range(n_part)]
    for i in range(n_part):
        x = random.random()*L
        y = random.random()*L
        vecx[i].append(x)
        vecy[i].append(y)
        arr.append((x,y))
    endo_rate = 10
    exo_rate = 10
    sigma_rate = 5
    freq = 100
    lamRecycleFreqRab4 = 20
    lamRecycleFreqTubular = 40
    recycleFreqTubular = 40 
    recycleFreqRab4 = 20
    rab4_efficiency = 0.4
    tubular_efficiency = 0.08
    img_names = []  
    for t in range(n_iter):
        if t%freq==0:
            fig, ax = plt.subplots()
            to_printx =[]
            to_printy =[]
            for i in range(n_part):
                to_printx.append(vecx[i][-1])
                to_printy.append(vecy[i][-1])
            ax.plot(to_printx, to_printy,'o')
            fig.savefig(f"images_my/{mechanism_type}/graph_{t}.png")
            plt.close(fig)
            img_names.append(f"images_my/{mechanism_type}/graph_{t}.png")
            nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
            while(nd>=n_part):
                nd = int(np.round(random.gauss(endo_rate,sigma_rate)))
            recycleFreqTubular = np.random.poisson(lam=lamRecycleFreqTubular)
            recycleFreqRab4 = np.random.poisson(lam=lamRecycleFreqRab4)
            recycled_particles_rab4 = Rab4Recycle(removed_particles_entering_rab4= nd*0.5, rab4_efficiency=Rab4Efficiency(n_part))
            recycled_particles_rab4_supposed = Rab4Recycle(removed_particles_entering_rab4= nd*0.5, rab4_efficiency=0.2)
            recycled_particles_tubular = TubularRecycle(removed_particles_entering_tubular=nd*0.5, tubular_efficiency=tubular_efficiency)
            recycled_particles_tubular_supposed = TubularRecycle(removed_particles_entering_tubular=nd*0.5, tubular_efficiency=0.4)
            removed_ind = set()
            while(len(removed_ind)<nd):
                x = random.randint(0, n_part-1)
                removed_ind.add(x)
            new_indices = []
            for i in range(n_part):
                if i not in removed_ind:
                    new_indices.append(i)
            tvecx=[]
            tvecy = []
            tarr=[]
            if len(new_indices) <= 50:
                rab4_efficiency = 0.5
            elif len(new_indices) >= 75:
                rab4_efficiency = 0.3
            for ind in new_indices:
                tvecx.append(vecx[ind])
                tvecy.append(vecy[ind])
                tarr.append(arr[ind])
            vecx = tvecx
            vecy = tvecy
            arr=tarr
            na = int(np.round(random.gauss(exo_rate,sigma_rate))) - (recycled_particles_rab4_supposed + recycled_particles_tubular_supposed)
            for i in range(na):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
            part_arr.append(n_part)
        elif (t+recycleFreqRab4)%freq==0:
            for i in range(recycled_particles_rab4):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
        elif (t+recycleFreqTubular)%freq==0:
            for i in range(recycled_particles_tubular):
                x = random.random()*L
                y = random.random()*L
                vecx.append([x])
                vecy.append([y])
                arr.append((x,y))
            n_part = len(vecx)
        if t%10 == 0:
            numbers.append(n_part)
    total_sum = 0
    for part in range(0, len(part_arr)):
        total_sum += part_arr[part]
    print(total_sum/len(part_arr))
    plt.figure(figsize=(10,10))
    plt.plot(numbers)
    plt.ylim(bottom=0)
    plt.title('Line Graph of E-cad molecules on surface vs Time')
    plt.xlabel('Time (Every 10th iteration)')
    plt.ylabel('Number of Particles')
    plt.savefig(f"graphs/{mechanism_type}/graph_{currentDate} {currentTime}.png")
    return img_names
def modelECAD():
    currentDate = datetime.now().date()
    currentTime = datetime.now().time()
    print("Generating random value for V between 0 and 1 ")
    n_iter = int(input("Enter the number of iterations :"))
    mechanism_type = str(input("Which Type of Mechanism do you want to run? A. Wild Type B. Rab11 C. RabX1 "))
    img_names = []
    if mechanism_type == "A" or mechanism_type == "a":
        mechanism_type = "wild_type"
        img_names = wildType(n_iter, mechanism_type)
    elif mechanism_type == "B" or mechanism_type == "b":
        mechanism_type = "rab11"
        img_names = Rab11(n_iter, mechanism_type)
    elif mechanism_type == "C" or mechanism_type == "c":
        mechanism_type = "rabx1"
        img_names = RabX1(n_iter, mechanism_type)
    print("pre-processing done!")
    codec = cv2.VideoWriter_fourcc(*"mp4v")
    out = cv2.VideoWriter(f"videos_my/{mechanism_type}/output {currentDate} {currentTime}.mp4", codec, 18, (640, 480))
    for i in range(len(img_names)):
        img = cv2.imread(img_names[i])
        out.write(img)
    plt.show()
    out.release()
    cv2.destroyAllWindows()
    print("done!")
# modelECAD()

# Plotting the inverse sigmoid function

# Define the inverse sigmoid function
def inverse_sigmoid(x):
    return 1 - 1 / (1 + np.exp(-x))

# Generate x values for plotting
x_values = np.linspace(-10, 10, 300)
y_values = inverse_sigmoid(x_values)

# Plotting
plt.figure(figsize=(10, 6))
plt.plot(x_values, y_values, label='Inverse Sigmoid Function')
plt.title('Rab4 Efficiency vs Number of Particles at Cell Surface')
plt.xlabel('Number of Particles at Cell Surface')
plt.ylabel('Rab4 Efficiency')
plt.grid(True)
# plt.gca().set_xlabel('')
# plt.gca().set_ylabel('')

# Remove x and y axis tick marks
plt.gca().set_xticks([])
plt.gca().set_yticks([])
plt.show()
