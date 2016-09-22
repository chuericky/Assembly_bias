### Stack the number densities of the galaxies with similar mass but different concentration quartiles.
### Updated: Aug 19, 2016

import numpy as np
import os, re, sys, glob
import matplotlib.pyplot as plt

quar_link = "/Users/rickyccy/Documents/Research_assemblybias/Data_multi/extension_work/galaxy/quartiles/"
result_link = "/Users/rickyccy/Documents/Research_assemblybias/result/Aug_16_2016/"

os.chdir(quar_link)

fd_num = 4      ### Number of fields
quar_num = 4    ### Number of quartiles

quar_gnum = [18806,18805,18805,18806]

#den_bin = ["3", "5", "7", "9", "11"]
#den_num = 5

for p in range(1):
    
    lens_pos = [0] * quar_num
    den_avg = [0] * quar_num
    
    print p

    for i in range(quar_num):      ### quar_num

        lens_pos_arr = [0] * fd_num
    
        for j in range(fd_num):
            lens_pos_arr[j] = np.loadtxt("W{}_iso_areabin{}_overdensity.dat".format(j+1,i+1), usecols=(1,2,3,4,5))
            #lens_pos_arr[j] = np.loadtxt("W{}_den{}Mpcbin{}_overdensity.dat".format(j+1,den_bin[p],i+1), usecols=(1,2,3,4,5))

        lens_pos[i] = np.vstack((lens_pos_arr[0], lens_pos_arr[1]))
        lens_pos[i] = np.vstack((lens_pos[i], lens_pos_arr[2]))
        lens_pos[i] = np.vstack((lens_pos[i], lens_pos_arr[3]))

        del lens_pos_arr

        den_avg[i] = np.average(lens_pos[i], axis = 0)


    d_avg = np.matrix(den_avg)

    de_avg = np.array(np.average(d_avg, axis = 0))[0]


    ### Jackknife
    jack_mat = [0] * quar_num
    sd_jack = [0] * quar_num

    ### Divide the galaxies into 50 jackknifes.
    jack_num = 50

    for i in range(quar_num):      ### quar_num
        m = lens_pos[i].T[0].size
        m_num = m/jack_num
        ### Number of jackknifes with 1 more galaxy.
        num_onemore = m % jack_num
        jack_arr = [0] * jack_num
    
        jack_arr[0] = m_num
        for k in range(1,jack_num - num_onemore):
            jack_arr[k] = jack_arr[k-1] + m_num
        for k in range(jack_num - num_onemore, jack_num):
            jack_arr[k] = jack_arr[k-1] + m_num + 1
        
        if (jack_arr[jack_num - 1] != m):
            print "Jackknife number isn't right!\n"
            exit(1)

        ### Sequence of jackknife
        jack_seq = np.random.choice(m,m,replace = False)

        if (np.where((np.bincount(jack_seq) != 1))[0].size != 0):
            print "Jackknife draw has repetition!\n"
            exit(1)

        ### Indices of galaxies to be eliminated in jackknife samples.
        jack_del = [0] * jack_num

        jack_del[0] = jack_seq[:jack_arr[0]]
        for k in range(1, jack_num):
            jack_del[k] = jack_seq[jack_arr[k-1]:jack_arr[k]]

        jack_B = [0] * jack_num
        jack_mat[i] = [0] * jack_num

        for k in range(jack_num):
            jack_B[k] = np.delete(lens_pos[i], jack_del[k], 0)
            jack_mat[i][k] = np.average(jack_B[k], axis = 0)

        jack_mat[i] = np.matrix(jack_mat[i])

        cov_jack = np.cov(jack_mat[i].T) * (jack_num - 1.)**2 / jack_num
        sd_jack[i] = np.power(np.diag(cov_jack),0.5)

        del jack_arr, jack_seq, jack_del, cov_jack


    ratjack_err = np.divide(sd_jack,de_avg)

    """
    a = np.arange(3,12,2)

    color = ["b", "g", "c", "r"]

    f, (ax1, ax2) = plt.subplots(1,2, figsize = (12,6))

    for i in range(quar_num):
        print ratjack_err[i]
        ax1.errorbar(a, den_avg[i], yerr=[sd_jack[i],sd_jack[i]], marker = "x", color = color[i])
        ax2.errorbar(a, np.divide(den_avg[i],de_avg), yerr=[ratjack_err[i],ratjack_err[i]], marker = "x", color = color[i])

    ax1.legend([r"Lowest $\delta$, 0% - 25%", "25% - 50%", "50% - 75%", "Highest $\delta$, 75% - 100%"], loc = 'best')
    ax1.set_xlim((2,12))
    ax2.set_xlim((2,12))
    ax1.set_xlabel(r"$r (h_{70}^{-1}\rm{Mpc})$")
    ax2.set_xlabel(r"$r (h_{70}^{-1}\rm{Mpc})$")
    ax1.set_ylabel(r"$N(r)$" + ", number of neighbors")
    ax2.set_ylabel(r"$\frac{N(r)}{\bar{N}}$")
    ax2.set_title(r"$\langle$" + "{}".format(den_bin[p]) + r"Mpc|\vec{p} \rangle$, Jackknife")
    ax2.hlines(xmin = 2,xmax=12,y=1,linestyle = '-.',color = 'k')
    ax2.set_ylim((0.9,1.1))
    plt.savefig(result_link + "Neighborcount_{}Mpc_jackknife.png".format(den_bin[p]))
    plt.close()
    """


    """
    for i in range(1):
        m = lens_pos[i].T[0].size
        B = [0] * m
        jack_mat[i] = [0] * m
    
        for j in range(m):
            B[j] = np.delete(lens_pos[i],j,0)
            jack_mat[i][j] = np.average(B[j],axis=0)
            #print i,j,jack_mat[i][j]

        jack_mat[i] = np.matrix(jack_mat[i])

        cov = np.cov(jack_mat[i].T) * m**2 / (m+1.)
        sd_mat = np.power(np.diag(cov),0.5)
        print sd_mat

        del B
    """


    ### Bootstrap
    rad_num = 5
    N_ran = 10000

    sd_navg = [0] * quar_num

    for i in range(quar_num):          ### quar_num
        print i
        N_size = lens_pos[i].T[0].size

        ### Bootstrap
        boot_num = [0] * N_ran
        n_avg = [0] * rad_num
    
        for j in range(rad_num):
            n_avg[j] = [0] * N_ran

        for j in range(N_ran):
            boot_num[j] = np.random.randint(0, N_size, size = N_size)

        for j in range(rad_num):
            for k in range(N_ran):
                n_avg[j][k] = np.average(np.array(lens_pos[i].T[j])[boot_num[k]])

        cov_navg = np.cov(n_avg)
        sd_navg[i] = np.power(np.diag(cov_navg),0.5)

        print cov_navg
        print sd_navg[i]

        del cov_navg

    rat_err = np.divide(sd_navg,de_avg)


    a = np.arange(3,12,2)

    color = ["b", "g", "c", "r"]
    plot_fac = [-0.4,-0.15,0.15,0.4]

    f, ((ax1, ax2, ax3), (ax4, ax5, ax6)) = plt.subplots(2,3, figsize = (20,12))

    for i in range(quar_num):
        print rat_err[i]
        den_y = np.divide(den_avg[i],de_avg)
        w_boot = np.power(np.divide(1,rat_err[i]),2)
        w_jack = np.power(np.divide(1,ratjack_err[i]),2)
        
        bfit_boot = np.polyfit(a, den_y, deg = 0, w = w_boot)
        bfit_jack = np.polyfit(a, den_y, deg = 0, w = w_jack)
        
        ax1.errorbar(a, den_avg[i], yerr=[sd_navg[i],sd_navg[i]], marker = "x", color = color[i])
        ax2.errorbar(a, den_y, yerr=[rat_err[i],rat_err[i]], marker = "x", color = color[i])
        ax3.errorbar(a + plot_fac[i], den_y - bfit_boot, yerr=[rat_err[i],rat_err[i]], marker = "x", color = color[i])
        ax4.errorbar(a, den_avg[i], yerr=[sd_jack[i],sd_jack[i]], marker = "x", color = color[i])
        ax5.errorbar(a, den_y, yerr=[ratjack_err[i],ratjack_err[i]], marker = "x", color = color[i])
        ax6.errorbar(a + plot_fac[i], den_y - bfit_jack, yerr=[ratjack_err[i],ratjack_err[i]], marker = "x", color = color[i])
        print i, i, np.divide(den_avg[i],de_avg)

    ax1.legend([r"Lowest area, 0% - 25%", "25% - 50%", "50% - 75%", "Highest area, 75% - 100%"], loc = 'best')
    #ax1.legend([r"Lowest $\delta$, 0% - 25%", "25% - 50%", "50% - 75%", r"Highest $\delta$, 75% - 100%"], loc = 'best')
    ax1.set_xlim((2,12))
    ax2.set_xlim((2,12))
    ax3.set_xlim((2,12))
    ax4.set_xlim((2,12))
    ax5.set_xlim((2,12))
    ax6.set_xlim((2,12))
    ax1.set_xlabel(r"$r (h_{70}^{-1}\rm{Mpc})$")
    ax2.set_xlabel(r"$r (h_{70}^{-1}\rm{Mpc})$")
    ax3.set_xlabel(r"$r (h_{70}^{-1}\rm{Mpc})$")
    ax1.set_ylabel(r"$N(r)$" + ", number of neighbors")
    ax2.set_ylabel(r"$\frac{N(r)}{\bar{N}}$")
    ax3.set_ylabel(r"$res(\frac{N(r)}{\bar{N}})$")
    ax4.set_xlabel(r"$r (h_{70}^{-1}\rm{Mpc})$")
    ax5.set_xlabel(r"$r (h_{70}^{-1}\rm{Mpc})$")
    ax6.set_xlabel(r"$r (h_{70}^{-1}\rm{Mpc})$")
    ax4.set_ylabel(r"$N(r)$" + ", number of neighbors")
    ax5.set_ylabel(r"$\frac{N(r)}{\bar{N}}$")
    ax6.set_ylabel(r"$res(\frac{N(r)}{\bar{N}})$")
    ax1.text(9,1000,"Bootstrap")
    ax4.text(9,1000,"Jackknife")
    ax2.text(9,0.905,"Bootstrap")
    ax5.text(9,0.905,"Jackknife")
    ax1.set_title("iso area")
    #ax2.set_title(r"$\langle \delta($" + "{}".format(den_bin[p]) + r"$Mpc)|\vec{p} \rangle$")
    ax2.hlines(xmin = 2,xmax=12,y=1,linestyle = '-.',color = 'k')
    ax5.hlines(xmin = 2,xmax=12,y=1,linestyle = '-.',color = 'k')
    ax3.hlines(xmin = 2,xmax=12,y=0,linestyle = '-.',color = 'k')
    ax6.hlines(xmin = 2,xmax=12,y=0,linestyle = '-.',color = 'k')
    ax2.set_ylim((0.8,1.2))
    ax3.set_ylim((-0.02,0.02))
    ax5.set_ylim((0.8,1.2))
    ax6.set_ylim((-0.02,0.02))
    ax1.set_ylim((0,35000))
    ax4.set_ylim((0,35000))
    plt.savefig(result_link + "Neighborcount_isoarea.png")
    plt.close()

    del den_avg, sd_navg, de_avg, sd_jack, lens_pos
