Program molecular_dynamics
    !! Molecular dynamics code for simulating a LJ cluster
    !! Output: file called 'msd_vs_time.out' containing: time, msd
    !! Output: file called 'temp_time.out' containing: time, temperature
    !! Output: file called 'energy.out' containing: time, energy
    !! Output: file called 'vmd.xyz' that is compatible with VMD to create movie
    !! Output: file called 'init_cond.out' containing initial state of particle in vmd format

    use mod_verlet
    implicit none
    real*8 tic, toc

    call cpu_time(tic)
    call setup_parameters
    call initial_condition
    call equilibrate
    call evolve

    call cpu_time(toc)
    write(*,*) "Total time taken by code = ",toc-tic

end program molecular_dynamics


