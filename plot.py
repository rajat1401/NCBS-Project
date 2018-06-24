import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.animation as animation


'''def bar_plot(values, dist):
    current_milli_time= lambda: int(round(time.time()*1000))

    while (int(round(time.time()*1000)< current_milli_time+2000):
        index= np.arange(len(dist))
        plt.bar (index, values)
        plt.xlabel('Genetic Codes', fontsize= 8)
        plt.ylabel('Relative Populations', fontsize= 8)
        #plt.xticks(index,   )
        plt.show()'''

fig= plt.figure()
plot= fig.add_subplot(1, 1, 1)

def smugsort(dist):
    a= sorted(range(len(dist)), key=lambda k: dist[k])
    return a


def bar_plot(i):
    data= open('sample', 'r').read()
    data= data.split('\n')
    values= []
    dist= []
    for each in data[:1]:
        dist= each.split(' ')

    print (dist)
    a= smugsort(dist)
    setgreen= 0
    isgreen= False
    dist.sort()
    for j in range(len(dist)):
        if(int(dist[j])== 0):#nice
            setgreen= j
            isgreen= True
            break

    color= ['blue']*len(dist)
    if(isgreen):
        color[setgreen]= 'green'

    index= np.arange(len(dist))
    data= data[i]
    values = []
    sortedvalues = []
    array= data.split(' ')
    for j in range(len(array)):
        values.append(float(array[j]))

    for j in range(len(array)):
        sortedvalues.append(values[a[j]])

    plot.clear()
    plot.bar(index, sortedvalues, color=color)
    '''for eachLine in data:
        if(len(eachLine)> 1):
            values= []
            sortedvalues= []
            array= eachLine.split(' ')
            for i in range(len(array)):
                values.append(float(array[i]))

            for i in range(len(array)):
                sortedvalues.append(values[a[i]])
            plot.clear()
            plot.bar(index, sortedvalues, color= color)'''


ani= animation.FuncAnimation(fig, bar_plot, interval= 500)#change this later!
plt.show()