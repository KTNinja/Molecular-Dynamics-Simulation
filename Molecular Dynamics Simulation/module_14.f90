Module mod_verlet
    implicit none

    integer npart
    real*8, allocatable :: pos(:,:),vel(:,:),acc(:,:),mass(:),mass_inv(:),pos_t0(:,:)
    real*8 pot, energy
    real*8 dt,total_time,curr_time, reset_temp_time, eq_time, final_time
    integer nsteps, nstep_reset, nstep_eq
    real*8 epsilon,sigma,sigma6,eps4,eps48, pi, cons
    real*8 density,temperature,box_len, temp_in
    real*8 E_avg,Esq_avg,E_std_dev
    real*8 V_C, sigma6_rc6, rc_2, msd
    integer flag_write, flag_vmd

    contains
    !##############################################

    ! Step 1: setting up parameters and allocating required arrays
    ! Reads from the input file md_LJ_fluid.inp
    subroutine setup_parameters
        implicit none
        real*8 mass_atom

        open(10,file="md_LJ_fluid.inp")
        read(10,*) npart
        read(10,*) dt
        read(10,*) reset_temp_time
        read(10,*) eq_time
        read(10,*) final_time
        read(10,*) mass_atom
        read(10,*) sigma
        read(10,*) epsilon
        read(10,*) density
        read(10,*) temperature
        read(10,*) flag_write
        read(10,*) flag_vmd
        close(10)

        allocate(pos(npart,3),vel(npart,3),acc(npart,3),mass(npart),mass_inv(npart), pos_t0(npart,3))
        mass=mass_atom
        mass_inv=1.d0/mass

        sigma6=sigma**6
        eps4=epsilon*4
        eps48=epsilon*48
        rc_2= 6.25*(sigma**2)
        sigma6_rc6 = sigma6/rc_2**3
        V_C = eps4*sigma6_rc6*(sigma6_rc6 - 1.d0)
        pi = 3.141592653
        cons = 2.d0/(3.d0*npart)

    end subroutine setup_parameters
    !##############################################

    !! Step 2: setting initial conditions
    !! initial position set to a simple cubic lattice
    !! initial velocities set to velocity corresponding to temperature of 0.71
    subroutine initial_condition
        implicit none
        integer i,j,k
        integer npart_side,part_num
        real*8 lattice_len

        npart_side=nint(1.d0*npart**(1/3.d0))   !! no. of particles along each direction
        box_len=(npart/density)**(1/3.d0)       !! length of simulation box
        lattice_len = box_len/real(npart_side)  !! length of the cubic unit cell

        write(*,*) "Initial conditions: placing particles in simple cubic lattice:"
        write(*,*) "total number of particles = ",npart
        write(*,*) "particles along each side = ",npart_side
        write(*,*) "box length = ",box_len
        write(*,*) "lattice length = ",lattice_len

        pos=0.d0
        !vel=0.d0
        call thermal_velocities(temperature)     !! initial velocity

        !! putting particles on a simple cubic lattice
        part_num=1
        do i=1,npart_side
            do j=1,npart_side
                do k=1,npart_side
                    pos(part_num,1)=(i-1)*lattice_len
                    pos(part_num,2)=(j-1)*lattice_len
                    pos(part_num,3)=(k-1)*lattice_len
                    part_num=part_num+1
                enddo
            enddo
        enddo



        !! Writing initial conditions in a vmd compatible xyz file
        open(13,file="init_cond.xyz")
        call write_vmd
        close(13)

        call compute_pot
        call compute_energy
        curr_time=0.0

        E_avg=0.d0
        Esq_avg=0.d0

    end subroutine initial_condition
    !##############################################

    !! Step 3: Equilibrating by calling velocity at every nstep_reset steps
    !! Writes vmd output if flag_vmd is set to 1 in the input file
    !! Writes temp and energy output if flag_write is set to 1 in the input file

    subroutine equilibrate
        implicit none
        integer i

        nstep_eq=nint(eq_time/dt)
        nstep_reset=nint(reset_temp_time/dt)

        open(18,file="temp_time.out")
        open(11,file="energy.out")
        open(13,file="vmd.xyz")

        do i=1,nstep_eq
            call compute_energy
            if(mod(i,nstep_reset)==1) call thermal_velocities(temperature)
            temp_in = cons*(energy-pot)
            if(flag_write==1 .and. mod(i,50)==1) call write_temp
            if(flag_write==1 .and. mod(i,50)==1) call write_energy
            if(flag_vmd==1 .and. mod(i,100)==1) call write_vmd
            call velocity_verlet
            call get_inside_the_box
            curr_time=curr_time+dt
        enddo

        close(18); close(11); close(13)

        pos_t0 = pos

    end subroutine equilibrate
    !##############################################

    !! Step 4: Evolution using velocity-Verlet
    !! Writes msd output if flag_write is set to 1 in the input file

    subroutine evolve
        implicit none
        integer i


        nsteps=nint(final_time/dt)
        curr_time = 0.d0
        open(10,file="msd_vs_time.out")

        do i=1,nsteps
            call compute_energy
            if(flag_write==1 .and. mod(i,5)==1) call write_output_msd
            call compute_msd
            call velocity_verlet
            call get_inside_the_box
            curr_time=curr_time+dt
        enddo
        !E_avg=E_avg/real(nsteps)
        !Esq_avg=Esq_avg/real(nsteps)
        !E_std_dev = dsqrt(Esq_avg-E_avg**2)
        !write(*,*) "dt (a.u.) standard deviation of E (a.u)",dt, E_std_dev
        close(10)

    end subroutine evolve
    !##############################################

    !! velocity verlet algorithm
    !! also calls for computing potential and acceleration
    subroutine velocity_verlet
        implicit none

        pos=pos+vel*dt+0.5*acc*dt*dt
        vel=vel+0.5*acc*dt
        call compute_pot
        vel=vel+0.5*acc*dt

    end subroutine velocity_verlet
    !##############################################

    !! calculates pairwise LJ potential and acceleration
    subroutine compute_pot
        implicit none
        integer i,j
        real*8 rijsq,vec(3),V_LJ,dVLJ_drij
        real*8 dpot_dx(3)

        pot=0.d0
        acc=0.d0
        do i=1,npart-1
            do j=i+1,npart
                call distance(i,j,rijsq,vec)
                if(rijsq<rc_2) then
                    call LJ_potential(rijsq,V_LJ,dVLJ_drij)
                    pot=pot+V_LJ
                    dpot_dx=dVLJ_drij*vec
                    acc(i,:)=acc(i,:)-mass_inv(i) * (dpot_dx)
                    acc(j,:)=acc(j,:)+mass_inv(j) * (dpot_dx)
                end if
            enddo
        enddo

    end subroutine compute_pot
    !##############################################

    !! calculates mean square displacement

    subroutine compute_msd
        implicit none

        integer i, k
        real*8 r_squared, vect(3)

        msd=0.d0
        do i=1,npart
            vect=pos(i,:)-pos_t0(i,:)
            do k=1,3
                if(vect(k) > 0.5*box_len) vect(k)= vect(k)-box_len
                if(vect(k) < -0.5*box_len) vect(k)= vect(k)+box_len
            end do
            r_squared=sum(vect*vect)
            msd=msd+r_squared
        end do
        msd=msd/npart

    end subroutine compute_msd
    !##############################################

    subroutine gaussian_random_number(rnd)

        implicit none

        double precision,intent(out)::rnd
        double precision U1,U2

        call random_number(U1)                              !!! Calculates a random number over a uniform interval [0,1]
        call random_number(U2)
        rnd=sqrt(-2*log(U1))*cos(2*pi*U2)

    end subroutine gaussian_random_number
    !#############################################

    !! uses the random number calculated above to calculate velocity

    subroutine thermal_velocities(temperature)

         implicit none
         double precision,intent(in) :: temperature
         double precision rnd,kb
         integer i,j

         kb=1.d0

         do i=1,npart
            do j=1,3
                call gaussian_random_number(rnd)            !! Random number with gaussian distribution with mean=0, sigma=1
                vel(i,j)=rnd*sqrt(kb*temperature/mass(i))   !! Random number with gaussian distribution with mean=0, sigma=sqrt(kb*temperature/mass(i))
            enddo
         enddo

    end subroutine thermal_velocities
    !##############################################

    !! calculates the distance squared between atoms i and j
    subroutine distance(i,j,rijsq,vec)
        !! input: i and j, that represent atoms i and j
        !! Output: vec(3) = pos(i,:)-pos(j,:) (vector from atom j to i)
        !! Output: rijsq: square of the distance between atoms i and j
        implicit none
        integer,intent(in) :: i,j
        real*8,intent(out) :: rijsq,vec(3)

        integer k

        do k=1,3
            vec(k)=pos(i,k)-pos(j,k)

            if(vec(k) > 0.5*box_len) vec(k)= vec(k)-box_len
            if(vec(k) < -0.5*box_len) vec(k)= vec(k)+box_len
        end do

        rijsq = (sum(vec*vec))

    end subroutine distance
    !##############################################

    subroutine get_inside_the_box

    !! Applies PBC such that if a particle leaves the box, its image entering the box is then considered.
    !! Box is considered cubic here with lengths from 0 to box_len along each direction
        implicit none
        integer i,j

        do i=1,npart                              !! Loop over all particles
           do j=1,3                               !! Loop along x,y,z direction

              if (pos(i,j)>box_len) then
                pos(i,j) = pos(i,j)-box_len
              else if (pos(i,j) <0) then
                pos(i,j) = pos(i,j)+box_len
              end if

           end do
        end do

    end subroutine get_inside_the_box
    !##############################################

    !! calculates LJ potential and 1/rij*dVLJ_dx
    !! Optimized
    subroutine LJ_potential(xsq,V_LJ,dVLJ_dx)
        !! Takes xsq=x_squares as the input
        !! Calculates V_LJ=LJ potential and dVLJ_dx=(1/x).derivative of potential
        implicit none
        real*8,intent(in)::xsq
        real*8,intent(out)::V_LJ,dVLJ_dx
        real*8 sig6_x6

        sig6_x6=sigma6/xsq**3

        V_LJ=(eps4*sig6_x6*(sig6_x6-1.d0)) - V_C
        dVLJ_dx=eps48/xsq*sig6_x6*(0.5-sig6_x6)

    end subroutine LJ_potential
    !##############################################

    !! writes mean square displacement for the current time step
    subroutine write_output_msd
        implicit none

        write(10,*) curr_time,msd

    end subroutine write_output_msd
    !##############################################

    !! writes energy for the current time step

    subroutine write_energy
        implicit none

        write(11,*) curr_time,energy

    end subroutine write_energy
    !##############################################

    !! writes temperature for the current time step

    subroutine write_temp
        implicit none

        write(18,*) curr_time, temp_in

    end subroutine write_temp
    !##############################################

    !! Writes vmd compatible file
    subroutine write_vmd
        implicit none
        integer i

        write(13,*) npart
        write(13,*)
        do i=1,npart
          write(13,*) "C",pos(i,:)
        enddo

    end subroutine write_vmd
    !##############################################

    !! calculates energy
    subroutine compute_energy
        implicit none
        integer i

        energy=pot
        do i=1,npart
          energy=energy+0.5*mass(i)*sum(vel(i,:)*vel(i,:))
        enddo

        !E_avg=E_avg+energy
        !Esq_avg=Esq_avg+energy**2

    end subroutine compute_energy
    !##############################################

end module mod_verlet
