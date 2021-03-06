from LIMBR import simulations, imputation, batch_fx

#20, 15, 15 missing and 5, 10, 15 NN

for i in range(1, 21):
    simulation = simulations.simulate(effect_size=4, p_miss=.4, lam_miss=15, amp_noise=1, rseed=i)
    simulation.generate_pool_map()
    simulation.write_output('twenty_miss_5_NN_'+str(i))
    simulation.write_output('twenty_miss_10_NN_'+str(i))
    simulation.write_output('twenty_miss_15_NN_'+str(i))
    simulation.write_output('thirty_miss_5_NN_'+str(i))
    simulation.write_output('thirty_miss_5_NN_'+str(i))
    simulation.write_output('thirty_miss_10_NN_'+str(i))
    simulation.write_output('thirty_miss_15_NN_'+str(i))
    simulation.write_output('forty_miss_5_NN_'+str(i))
    simulation.write_output('forty_miss_10_NN_'+str(i))
    simulation.write_output('forty_miss_15_NN_'+str(i))
