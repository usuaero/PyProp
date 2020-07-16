module prop_m
    use myjson_m
    use section_m
    implicit none
    
    type motor_t
        integer :: ID
        character(20) :: name
        real :: kv  ! kv rpm/V
        real :: Io ! No Load Current (Amps)
        real :: resistance ! Resistance (Ohms)
        real :: gearratio !gear ratio
        real :: rpm
        real :: torque !Torque (N-m)
        real :: current !current (Amps)
        real :: voltage !voltage (Volts)
        real :: brakepower !Break Power (Watts)
    end type motor_t

    type battery_t
        integer :: ID
        character(20) :: name
        integer :: ncells
        real :: capacity !mAh
        real :: Eo     !No Load Voltage
        real :: resistance !Resistance (Ohms)
        real :: crating !C-rating
        real :: voltage  !voltage (Volts)
        real :: current  !current (Amps)
        real :: power    ! power (Watts)
        real :: endurance !(minutes)
    end type battery_t


    type speedcontrol_t
        integer :: ID
        character(20) :: name
        real :: resistance !Ohms
        real :: eta !efficiency
    end type speedcontrol_t
    
    type prop_t
        integer :: ID
        character(10) :: name
        character(100) :: master_filename

        type(json_file) :: json    !the JSON structure read from the file:

        integer :: verbose !=1 for write info, =0 for don't write info
        
        integer :: nblades
        real :: diameter
        real :: pitch
        real :: hubradius
        real :: chord_1
        real :: chord_2
        character(3) :: rotation
        character(5) :: side

        integer :: nSec
        type(section_t),pointer,dimension(:) :: sec
        
        character(100) :: af1_text,af2_text !original text from file
        type(airfoil_t),pointer :: af1
        type(airfoil_t),pointer :: af2


        type(motor_t) :: motor
        type(battery_t) :: battery
        type(speedcontrol_t) :: esc
        
        real :: velocity
        real :: density
        
        
        real :: rpm
        real :: omega !rad/s
        real :: Thrust,Torque,Power
        real :: C_thrust,C_torque,C_power
        
        !possible files
        character(100) :: f_pitch, f_chord, f_root_airfoil, f_tip_airfoil

        !Other Geometry Specs
        real :: root(3)
        real :: tip(3)
        real :: span
        real :: area

        real :: throttle
        real :: eta !system efficiency

    end type prop_t
    
        character(20) :: run_type !needed to know how to calculate sytem properties in prop_forces

        !Airfoil Database
        character(200) :: DB_Airfoil
        type(airfoil_t),pointer,dimension(:) :: airfoils

contains

!-----------------------------------------------------------------------------------------------------------
subroutine prop_allocate(t)
    type(prop_t) :: t
!    if(size(t%sec).gt.1) call prop_deallocate(t)

    allocate(t%sec(t%nSec))
end subroutine prop_allocate

!-----------------------------------------------------------------------------------------------------------
subroutine prop_deallocate(t)
    type(prop_t) :: t
    deallocate(t%sec)
    call t%json%destroy()
end subroutine prop_deallocate

!-----------------------------------------------------------------------------------------------------------
subroutine prop_set_defaults(t)
    type(prop_t) :: t    

    t%rpm = 10000.0

end subroutine prop_set_defaults

!-----------------------------------------------------------------------------------------------------------
subroutine prop_load_json(t,error)
    type(prop_t) :: t
    type(json_value),pointer :: j_afprop
    character(len=:),allocatable :: cval
    integer :: error,iairfoil,loc

    call json_initialize()

!    write(*,*) 'reading input file: ',t%master_filename
    call t%json%load_file(filename = t%master_filename); call json_check()
    loc = index(t%master_filename,'.json')
    t%master_filename = t%master_filename(1:loc-1) !deletes the .json file extension

    write(*,*) 'Reading Input File...'

!    call t%json%print_file()

    call t%json%get('airfoil_DB',                        cval);     call json_check(); DB_Airfoil = trim(cval)
    write(*,*) 'Airfoil database located at: ',trim(DB_Airfoil)

    t%velocity = json_file_required_real(t%json,'condition.velocity')
    t%density = json_file_optional_real(t%json,'condition.density',1.225)

    call t%json%get('propeller.name',             cval);          call json_check(); t%name = trim(cval)
    call t%json%get('propeller.#blades',    t%nblades);     call json_check()
    call t%json%get('propeller.diameter',   t%diameter);    call json_check()
    call t%json%get('propeller.pitch',      t%pitch);       call json_check()
    call t%json%get('propeller.hub_radius', t%hubradius);   call json_check()
    call t%json%get('propeller.root_chord', t%chord_1);     call json_check()
    call t%json%get('propeller.tip_chord',  t%chord_2);     call json_check()
    call t%json%get('propeller.rotation',         cval);          call json_check(); t%rotation = trim(cval)
    call t%json%get('propeller.grid',       t%nSec);        call json_check()

    allocate(airfoils(2))
    iairfoil = 0

    !Root Airfoil
    iairfoil = iairfoil + 1
    call t%json%get('propeller.root_airfoil.name',                    cval); call json_check()
    t%af1_text = trim(cval)
    airfoils(iairfoil)%name = trim(t%af1_text)
    call t%json%get('propeller.root_airfoil.properties',  j_afprop);
    if(json_failed()) then !Read from airfoil database
        call json_clear_exceptions()
        call prop_load_airfoil(t,iairfoil,'',0)
    else !Read from local file
        call prop_load_airfoil(t,iairfoil,'propeller.root_airfoil.',1)
    end if
    t%af1 => airfoils(iairfoil)

    !Tip Airfoil
    iairfoil = iairfoil + 1
    call t%json%get('propeller.tip_airfoil.name',                    cval); call json_check()
    t%af2_text = trim(cval)
    airfoils(iairfoil)%name = trim(t%af2_text)
    call t%json%get('propeller.tip_airfoil.properties',  j_afprop);
    if(json_failed()) then !Read from airfoil database
        call json_clear_exceptions()
        call prop_load_airfoil(t,iairfoil,'',0)
    else !Read from local file
        call prop_load_airfoil(t,iairfoil,'propeller.tip_airfoil.',1)
    end if
    t%af2 => airfoils(iairfoil)


    t%f_pitch         = 'none'
    t%f_chord         = 'none'
    t%f_root_airfoil  = 'none'
    t%f_tip_airfoil   = 'none'
    t%side = 'right'
    if(t%rotation.eq.'CW') t%side = 'left'

    call t%json%get('motor.name',                 cval);                  call json_check(); t%motor%name = trim(cval)
    call t%json%get('motor.kv',             t%motor%kv);            call json_check()
    call t%json%get('motor.no_load_current',t%motor%Io);            call json_check()
    call t%json%get('motor.resistance',     t%motor%resistance);    call json_check()
    call t%json%get('motor.gear_ratio',     t%motor%gearratio);     call json_check()

    call t%json%get('battery.name',                      cval);                  call json_check(); t%battery%name = trim(cval)
    call t%json%get('battery.#cells',              t%battery%ncells);      call json_check()
    call t%json%get('battery.cell_capacity',       t%battery%capacity);    call json_check()
    call t%json%get('battery.cell_no_load_voltage',t%battery%Eo);          call json_check()
    call t%json%get('battery.cell_impedance',      t%battery%resistance);  call json_check()
    call t%json%get('battery.c_rating',            t%battery%crating);     call json_check()

    call t%json%get('esc.name',              cval);               call json_check(); t%esc%name = trim(cval)
    call t%json%get('esc.resistance',  t%esc%resistance);   call json_check()

end subroutine prop_load_json

!-----------------------------------------------------------------------------------------------------------
subroutine prop_load_airfoil(t,i,prefix,local)
    type(prop_t) :: t
    type(json_file) :: f_json    !the JSON structure read from the file:
    character(len=*) :: prefix
    integer :: local
    integer :: i,ios
    character(100) :: fn,datafilename
    character(len=:),allocatable :: cval
    
    if(local.eq.1) then
        fn = trim(t%master_filename)//'.json'
    else
        fn = trim(adjustl(DB_Airfoil))//'/'//trim(adjustl(airfoils(i)%name))//'.json'
    end if

    write(*,*) 'Reading airfoil properties from file: ',trim(fn)
    call f_json%load_file(filename = trim(fn)); call json_check()
!    call f_json%print_file()

    call f_json%get(trim(prefix)//'properties.type', cval);  call json_check(); airfoils(i)%properties_type = trim(cval)

    select case (airfoils(i)%properties_type)
        case ('linear')
            airfoils(i)%aL0 = json_file_required_real(f_json,trim(prefix)//'properties.alpha_L0');
            airfoils(i)%CLa = json_file_required_real(f_json,trim(prefix)//'properties.CL_alpha'); 
            airfoils(i)%CmL0 = json_file_required_real(f_json,trim(prefix)//'properties.Cm_L0');
            airfoils(i)%Cma = json_file_required_real(f_json,trim(prefix)//'properties.Cm_alpha'); 
            airfoils(i)%CD0 = json_file_required_real(f_json,trim(prefix)//'properties.CD0');
            airfoils(i)%CD0L = json_file_required_real(f_json,trim(prefix)//'properties.CD0_L');
            airfoils(i)%CD0L2 = json_file_required_real(f_json,trim(prefix)//'properties.CD0_L2');
            airfoils(i)%CLmax = json_file_optional_real(f_json,trim(prefix)//'properties.CL_max',-1.0);
            airfoils(i)%has_data_file = 0
        case ('datafile')
            call f_json%get(trim(prefix)//'properties.filename', cval); call json_check(); datafilename = trim(cval)
            datafilename = trim(adjustl(DB_Airfoil))//'/'//trim(adjustl(datafilename))
            call af_create_from_data_file(airfoils(i),datafilename)
            airfoils(i)%has_data_file = 1
        case default
            write(*,*) 'Invalid airfoil properties type! Aborting program.'
            stop
    end select

    write(*,*) 'Loaded airfoil: ',airfoils(i)%name
end subroutine prop_load_airfoil

!-----------------------------------------------------------------------------------------------------------
subroutine prop_write_json_file(t,filename)
    implicit none
    type(prop_t) :: t
    type(json_value),pointer    :: p_root, p_condition, p_prop, p_motor, p_batt, p_esc, p_run
    character(100) :: filename
    integer :: ios,iairfoil,i,iunit

    iairfoil = 0
    ios = 0

    !root:
    call json_value_create(p_root)           ! create the value and associate the pointer
    call to_object(p_root,trim(filename))    ! add the file name as the name of the overall structure
    call json_value_add(p_root, 'case', 'ABCDEFG')

    write(*,'(A)') ''
    write(*,'(A)') 'initialize the structure...'

    !condition structure:
    call json_value_create(p_condition)             !an object
    call to_object(p_condition,'condition')
    call json_value_add(p_root, p_condition)
    call add_real_to_json_obj(  p_condition, 'velocity', t%velocity, 'm/s')
    call add_real_to_json_obj(  p_condition, 'density',  t%density,  'kg/m^3')
    nullify(p_condition)

    !propeller structure:
    call json_value_create(     p_prop)
    call to_object(             p_prop, 'propeller')
    call json_value_add(p_root, p_prop)
    call json_value_add(        p_prop, 'name',         trim(t%name))
    call add_int_to_json_obj(   p_prop, '#blades',      t%nblades,   '')
    call add_real_to_json_obj(  p_prop, 'diameter',     t%diameter,  'm')
    call add_real_to_json_obj(  p_prop, 'pitch',        t%pitch,     'm')
    call add_real_to_json_obj(  p_prop, 'hub_radius',   t%hubradius, 'm')
    call add_real_to_json_obj(  p_prop, 'root_chord',   t%chord_1,   'm')
    call add_real_to_json_obj(  p_prop, 'tip_chord',    t%chord_2,   'm')
    call json_value_add(        p_prop, 'rotation',     t%rotation)
    call json_value_add(        p_prop, 'root_airfoil', t%af1_text)
    call json_value_add(        p_prop, 'tip_airfoil',  t%af2_text)
    call add_int_to_json_obj(   p_prop, 'grid',         t%nSec,       '')
    nullify(p_prop)

    !motor structure:
    call json_value_create(     p_motor)
    call to_object(             p_motor,'motor')
    call json_value_add(p_root, p_motor)
    call json_value_add(        p_motor, 'name',            trim(t%motor%name))
    call add_real_to_json_obj(  p_motor, 'kv',              t%motor%kv,        'rpm/v')
    call add_real_to_json_obj(  p_motor, 'no_load_current', t%motor%Io,        'A')
    call add_real_to_json_obj(  p_motor, 'resistance',      t%motor%resistance,'Ohms')
    call add_real_to_json_obj(  p_motor, 'gear_ratio',      t%motor%gearratio, '')
    nullify(p_motor)

    !battery structure:
    call json_value_create(     p_batt)
    call to_object(             p_batt,'battery')
    call json_value_add(p_root, p_batt)
    call json_value_add(        p_batt, 'name',                 trim(t%battery%name))
    call add_int_to_json_obj(   p_batt, '#cells',               t%battery%ncells,    '')
    call add_real_to_json_obj(  p_batt, 'cell_capacity',        t%battery%capacity,  'mAh')
    call add_real_to_json_obj(  p_batt, 'cell_no_load_voltage', t%battery%Eo,        'V')
    call add_real_to_json_obj(  p_batt, 'cell_impedance',       t%battery%resistance,'Ohms')
    call add_real_to_json_obj(  p_batt, 'c_rating',             t%battery%crating,   '')
    nullify(p_batt)

    !ESC structure:
    call json_value_create(     p_esc)
    call to_object(             p_esc,'esc')
    call json_value_add(p_root, p_esc)
    call json_value_add(        p_esc, 'name',          trim(t%esc%name))
    call add_real_to_json_obj(  p_esc, 'resistance',    t%esc%resistance,'Ohms')
    nullify(p_esc)

    write(*,'(A)') 'writing file '//trim(filename)//'...'
    open(newunit=iunit, file=filename, status='REPLACE')
    call json_print(p_root,iunit)
    close(iunit)
    call json_destroy(p_root)
    write(*,'(A)') 'done.'

end subroutine prop_write_json_file

!-----------------------------------------------------------------------------------------------------------
subroutine add_real_to_json_obj(me, name, value, units)
    implicit none
    type(json_value),pointer :: me, var
    character(len=*),intent(in) :: name, units
    real, intent(in) :: value
    nullify(var)
    call json_value_create(var)
    call to_object(var,name)
    call json_value_add(var, 'value', value)
    call json_value_add(var, 'units', trim(units))
    call json_value_add(me, var)
    nullify(var)
end subroutine add_real_to_json_obj
!-----------------------------------------------------------------------------------------------------------
subroutine add_int_to_json_obj(me, name, value, units)
    implicit none
    type(json_value),pointer :: me, var
    character(len=*),intent(in) :: name, units
    integer, intent(in) :: value
    nullify(var)
    call json_value_create(var)
    call to_object(var,name)
    call json_value_add(var, 'value', value)
    call json_value_add(var, 'units', trim(units))
    call json_value_add(me, var)
    nullify(var)
end subroutine add_int_to_json_obj

!-----------------------------------------------------------------------------------------------------------
subroutine prop_init_setup(t)
    type(prop_t) :: t
    integer :: isec

    write(*,*) 'Setting up case...'

    call prop_allocate(t)
    call prop_setup(t)
    call prop_write_attributes(t)

end subroutine prop_init_setup

!-----------------------------------------------------------------------------------------------------------
subroutine prop_run_rpm(t,rpm)
    type(prop_t) :: t
    real :: rpm
    run_type = 'rpm'
    t%rpm = rpm
    write(*,*) 'Target RPM = ',t%rpm

    call prop_solve(t,2)

end subroutine prop_run_rpm

!-----------------------------------------------------------------------------------------------------------
subroutine prop_run_throttle(t,throttle)
    type(prop_t) :: t
    integer :: i
    real :: throttle,rpm_error
    
    run_type = 'throttle'
    t%throttle = throttle
    rpm_error = 0.1
!    t%rpm = 10000.0
    t%motor%rpm = 0.0
    write(*,*) 'Target throttle = ',t%throttle

!    do i=1,10
    write(*,*)
    write(*,*) '  Propeller_RPM             Motor_RPM                 Error'
    do while (abs(t%rpm-t%motor%rpm*t%motor%gearratio) > rpm_error)
        write(*,*) t%rpm,t%motor%rpm*t%motor%gearratio,abs(t%rpm-t%motor%rpm*t%motor%gearratio)
        call prop_solve(t,0)
        t%rpm = t%rpm - 0.3*(t%rpm - t%motor%rpm*t%motor%gearratio)
    end do
    call prop_solve(t,2)
end subroutine prop_run_throttle

!-----------------------------------------------------------------------------------------------------------
subroutine prop_run_thrust(t,thrust)
    type(prop_t) :: t
    integer :: i
    real :: thrust,thrust_error,C0
    
    run_type = 'thrust'
    thrust_error = 0.001 !0.01=1.0%
!    t%rpm = 10000.0
    t%Thrust = 0.0
    write(*,*) 'Target thrust = ',thrust

!    do i=1,10
    write(*,*)
    write(*,*) '  RPM                       Thrust [N]'
    do while (abs(t%Thrust-thrust) > thrust_error)
        call prop_solve(t,0)
        write(*,*) t%rpm,t%Thrust
!        C0 = t%Thrust/(t%rpm)**2 !assumes quadratic and passes through origin
!        t%rpm = sqrt(thrust/C0)
        t%rpm = t%rpm + 100.0*(thrust-t%Thrust) !relaxation
end do
    call prop_solve(t,2)
end subroutine prop_run_thrust

!-----------------------------------------------------------------------------------------------------------
subroutine prop_solve(t,iprint)
    type(prop_t) :: t
    type(section_t),pointer :: si
    integer :: isec,iter,iprint
    real :: ei,ei_error,x0,x1,guess
    real :: time1,time2
!    call cpu_time(time1)


    ei_error = 0.000001
    guess = 0.1
    t%omega = t%rpm*2.0*pi/60.0

    do isec=t%nSec,1,-1 !backwards is more stable because each case starts with initial guess of previous.
        si => t%sec(isec)
        si%einf = atan(t%velocity/(t%omega*si%rc))

        ei = guess
        x0 = guess
        x1 = 1.1*guess
        iter = 0
!        write(*,*) isec,t%nblades, t%diameter, t%sec(t%nSec)%beta,si%rc,si%chord_c,si%einf,si%beta,si%aL0,prop_function(t,si,ei)

        !Use the secant method to solve for ei
        do while (abs(prop_function(t,si,ei)).gt.ei_error)
            iter = iter + 1
            ei = x1 - prop_function(t,si,x1)*(x1 - x0)/(prop_function(t,si,x1) - prop_function(t,si,x0))
            x0 = x1
            x1 = ei
!            write(*,*) iter,ei,prop_function(t,si,ei)
        end do
        si%ei = ei
        guess = ei
!        write(*,*) 'iterations: ',iter
    end do
    
    call prop_forces(t,iprint)
    
!    call cpu_time(time2)
!    write(*,*) 'CPU time to solve (sec): ',time2-time1
end subroutine prop_solve

!-----------------------------------------------------------------------------------------------------------
real function prop_function(t,si,ei)
    type(prop_t) :: t
    type(section_t),pointer :: si
    real :: ei
    integer :: k
    real :: dp,betat,r,cb,einf,beta,aL0
    
    k = t%nblades
    dp = t%diameter
    betat = t%sec(t%nSec)%beta

    r = si%rc
    cb = si%chord_c
    einf = si%einf
    beta = si%beta
    aL0 = si%aL0
    si%alpha = beta - einf - ei + aL0
    
    prop_function = real(k)*cb*sec_CL(si)/16.0/r - acos(exp(-k*(1.0-2.0*r/dp)/2.0/sin(betat)))*tan(ei)*sin(einf+ei)

end function prop_function

!-----------------------------------------------------------------------------------------------------------
subroutine prop_forces(t,iprint)
    type(prop_t) :: t
    type(section_t),pointer :: si
    integer :: isec,iprint
    real :: ei,einf,r,cb,dr,denom
    120 Format(8ES25.13)

    t%Thrust = 0.0
    t%Torque = 0.0
    t%Power = 0.0

    do isec=1,t%nSec
        si => t%sec(isec)
        ei = si%ei
        einf = si%einf
        r = si%rc
        cb = si%chord_c
        dr = abs(si%r2-si%r1) !left side is backwards
        !alpha is already set from the solver

        t%Thrust = t%Thrust + cb*  (r*cos(ei)/cos(einf))**2*(sec_CL(si)*cos(einf+ei) - sec_CD(si)*sin(einf+ei))*dr
        t%Torque = t%Torque + cb*r*(r*cos(ei)/cos(einf))**2*(sec_CD(si)*cos(einf+ei) + sec_CL(si)*sin(einf+ei))*dr
    end do
    
    t%Thrust = t%Thrust*t%nblades*t%density*t%omega**2/2.0
    t%Torque = t%Torque*t%nblades*t%density*t%omega**2/2.0
    t%Power = t%Torque*t%omega
    
    denom = t%density*(t%omega/2.0/pi)**2*t%diameter**4

    t%C_thrust = t%Thrust/denom
    t%C_torque = t%Torque/denom/t%diameter
    t%C_power = t%Power/denom/(t%omega/2.0/pi)/t%diameter

    if(iprint.gt.0) then
        !write forces to the screen
        write(*,*) '------------ Propeller Solution ---------------'
        write(*,*) '     Thrust(N)                Torque(N-m)              Power(W)                 C_Thrust                 &
                        &C_Torque                 C_Power'
        write(*,120) t%Thrust,t%Torque,t%Power,t%C_thrust,t%C_torque,t%C_power
    end if

    !update motor/battery/prop combo
    t%motor%torque      = t%Torque/t%motor%gearratio
    t%motor%current     = pi/30.0*t%motor%kv*t%motor%torque + t%motor%Io
    t%battery%voltage   = t%battery%Eo * real(t%battery%ncells) - t%motor%current*t%battery%resistance * real(t%battery%ncells)

    if(run_type.eq.'throttle') then
        t%battery%current   = t%throttle*t%motor%current
        t%esc%eta           = 1.0-0.078*(1.0-t%throttle)
        t%motor%voltage     = t%esc%eta*t%throttle*t%battery%voltage - t%motor%current*t%esc%resistance
        t%motor%rpm         = t%motor%kv*(t%motor%voltage - t%motor%current*t%motor%resistance)
    else
        t%motor%rpm = t%rpm/t%motor%gearratio
        t%motor%voltage = t%motor%rpm/t%motor%kv + t%motor%current*t%motor%resistance
        t%throttle = (t%motor%voltage + t%motor%current*t%esc%resistance)/t%battery%voltage
        t%esc%eta = 1.0 !Fix This
        t%battery%current   = t%throttle*t%motor%current
    end if
    t%motor%brakepower  = pi/30.0*t%motor%torque*t%motor%rpm
    t%battery%endurance = t%battery%capacity*real(t%battery%ncells)/1000.0/t%battery%current*60.0
    t%battery%power = t%battery%current * t%battery%voltage
    t%eta = t%motor%brakepower/t%battery%power
    

    
    if(iprint.gt.1) then
        write(*,*) '------------ System Properties --------------'
        write(*,*) '      Motor Torque (Nm) = ',t%motor%torque
        write(*,*) '      Motor Current (A) = ',t%motor%current
        write(*,*) '      Motor Voltage (V) = ',t%motor%voltage
        write(*,*) '              Motor RPM = ',t%motor%rpm
        write(*,*) '  Motor Brake Power (W) = ',t%motor%brakepower
        write(*,*)
        write(*,*) '    Battery Current (A) = ',t%battery%current
        write(*,*) '    Battery Voltage (V) = ',t%battery%voltage
        write(*,*) '      Battery Power (W) = ',t%battery%power
        write(*,*) 'Battery Endurance (min) = ',t%battery%endurance
        write(*,*)
        write(*,*) '      System Efficiency = ',t%eta
        write(*,*) '               Throttle = ',t%throttle
        write(*,*) '             Thrust (N) = ',t%Thrust
        write(*,*)
    end if

end subroutine prop_forces

!-----------------------------------------------------------------------------------------------------------
subroutine prop_distributions(t)
    type(prop_t) :: t
    character(100) :: filename
    integer :: isec,ierror
    type(section_t),pointer :: si
    real :: ei, einf, r, cb, dr, my_thrust, my_torque
    120 Format(20ES25.13)

    filename = trim(adjustl(t%master_filename))//'_distributions.txt'
    open(unit = 10, File = filename, action = 'write', iostat = ierror)
    write(10,*) 'ControlPoint(x)          ControlPoint(y)          ControlPoint(z)          Chord(m)                 &
                &Twist(deg)               alpha(deg)               CL                       CD                       &
                &Cm                       area                     Thrust(N)                Torque(Nm)'
    do isec=1,t%nSec
        si => t%sec(isec)
        ei = si%ei
        einf = si%einf
        r = si%rc
        cb = si%chord_c
        dr = abs(si%r2-si%r1) !left side is backwards
        !alpha is already set from the solver

        my_thrust = cb*  (r*cos(ei)/cos(einf))**2*(sec_CL(si)*cos(einf+ei) - sec_CD(si)*sin(einf+ei))*dr
        my_torque = cb*r*(r*cos(ei)/cos(einf))**2*(sec_CD(si)*cos(einf+ei) + sec_CL(si)*sin(einf+ei))*dr
        
        my_thrust = my_thrust*t%density*t%omega**2/2.0
        my_torque = my_torque*t%density*t%omega**2/2.0
    
        write(10,120) si%PC(:),si%chord_c,si%twist*180.0/pi,si%alpha*180.0/pi,sec_CL(si),sec_CD(si),sec_Cm(si),si%ds,&
                    & my_thrust,my_torque
    end do
    close(10)
end subroutine prop_distributions

!-----------------------------------------------------------------------------------------------------------
subroutine prop_setup(t)
    type(prop_t) :: t
    integer :: isec
    type(section_t),pointer :: si
    real :: r1,r2,rc,aL0,lambda,beta,temp
    real :: qvec(3),nvec(3),avec(3),start(3),dtheta,percent_1,percent_2,percent_c,chord_1,chord_2,span
    real :: my_twist, twist1, twist2
    type(dataset_t) :: data_pitch,data_chord

    !allocate arrays from files
    if(t%f_chord .ne. 'none') then
        call ds_create_from_file(data_chord,t%f_chord,2)
    end if
    if(t%f_pitch .ne. 'none') then
        call ds_create_from_file(data_pitch,t%f_pitch,2)
    end if

    t%span = t%diameter/2.0 - t%hubradius
    t%root = 0.0
    t%root(2) = t%hubradius
    if(t%side.eq.'left') t%root(2) = -t%hubradius
    start = t%root

    dtheta = pi/real(t%nSec)
t%area = 0.0
span = 0.0

    do isec=1,t%nSec

        si => t%sec(isec)
        
        percent_1 = 0.5*(1.0-cos(dtheta*real(isec-1)))
        percent_2 = 0.5*(1.0-cos(dtheta*real(isec)))
        percent_c = 0.5*(1.0-cos(dtheta*(real(isec)-0.5)))
        if(t%side.eq.'left') then !must handle differently for opposite motion
            temp = percent_1
            percent_1 = percent_2
            percent_2 = temp
        end if
        
        r1 = t%hubradius + percent_1*(t%span)
        r2 = t%hubradius + percent_2*(t%span)
        rc = t%hubradius + percent_c*(t%span)

        !data from wing file
        if(t%chord_2 >= 0.0) then !linear
            chord_1 = t%chord_1 + percent_1*(t%chord_2 - t%chord_1)
            chord_2 = t%chord_1 + percent_2*(t%chord_2 - t%chord_1)
        else !elliptic wing
            chord_1 = t%diameter*0.075*sqrt(1.0-(2.0*r1/t%diameter)**2)
            chord_2 = t%diameter*0.075*sqrt(1.0-(2.0*r2/t%diameter)**2)
        end if
        
        aL0 = t%af1%aL0 + percent_c*(t%af2%aL0 - t%af1%aL0)
        lambda = 2.0*pi*rc*(t%pitch - 2.0*pi*rc*tan(aL0))/(2.0*pi*rc + t%pitch*tan(aL0)) !aerodynamic pitch (length)
        si%beta = atan(lambda/(2.0*pi*rc)) !aerodynamic pitch angle
        my_twist = si%beta + aL0 !geometric pitch angle
        si%aL0 = aL0
        
        !For geometry purposes
        aL0 = t%af1%aL0 + percent_1*(t%af2%aL0 - t%af1%aL0)
        lambda = 2.0*pi*r1*(t%pitch - 2.0*pi*r1*tan(aL0))/(2.0*pi*r1 + t%pitch*tan(aL0))
        beta = atan(lambda/(2.0*pi*r1))
        twist1 = beta + aL0

        aL0 = t%af1%aL0 + percent_2*(t%af2%aL0 - t%af1%aL0)
        lambda = 2.0*pi*r2*(t%pitch - 2.0*pi*r2*tan(aL0))/(2.0*pi*r2 + t%pitch*tan(aL0))
        beta = atan(lambda/(2.0*pi*r2))
        twist2 = beta + aL0

        si%dihedral1 = 0.0
        si%dihedral2 = 0.0

        ! data from input files
        if(t%f_chord .ne. 'none') then
            chord_1 = t%chord_1 * ds_linear_interpolate_col(data_chord,percent_1,1,2)
            chord_2 = t%chord_1 * ds_linear_interpolate_col(data_chord,percent_2,1,2)
        end if
        if(t%f_pitch .ne. 'none') then !fix this. Need to update si%beta for twist file input
            my_twist   = ds_linear_interpolate_col(data_pitch,percent_c,1,2)*pi/180.0
            twist1     = ds_linear_interpolate_col(data_pitch,percent_1,1,2)*pi/180.0
            twist2     = ds_linear_interpolate_col(data_pitch,percent_2,1,2)*pi/180.0
            aL0 = t%af1%aL0 + percent_c*(t%af2%aL0 - t%af1%aL0)
            si%beta = my_twist - aL0 !aerodynamic pitch angle
        end if
        
        !for geometry purposes
        si%twist = my_twist
        si%twist1 = twist1
        si%twist2 = twist2
        
        !operate!
        qvec(1) = 0.0; qvec(2) = 1.0; qvec(3) = 0.0
        nvec(1) = 0.0; nvec(2) = 0.0; nvec(3) =-1.0
        avec(1) =-1.0; avec(2) = 0.0; avec(3) = 0.0

        call math_rot_y(nvec,my_twist)
        call math_rot_y(avec,my_twist)
        if(t%side.eq.'left') then !must handle differently for opposite motion
            qvec(2) = -qvec(2)
            nvec(2) = -nvec(2)
            avec(2) = -avec(2)
        end if

        si%r1 = r1
        si%r2 = r2
        si%rc = rc
        si%percent_1 = percent_1
        si%percent_2 = percent_2
        si%percent_c = percent_c
        si%chord_1 = chord_1
        si%chord_2 = chord_2
        si%chord_c = 2.0/3.0*(chord_1**2 + chord_1*chord_2 + chord_2**2)/(chord_1 + chord_2)
        si%ds = abs(0.5*(chord_1+chord_2)*(percent_2-percent_1)*t%span)
        si%un(:) = nvec(:)
        si%ua(:) = avec(:)
        call math_cross_product(avec(:),nvec(:),si%us)
        
        if(t%side.eq.'left') then !must handle differently for opposite motion
            si%P2(:) = start(:)
            si%P1(:) = start(:) + qvec(:)*(percent_1-percent_2)*t%span
        else
            si%P1(:) = start(:)
            si%P2(:) = start(:) + qvec(:)*(percent_2-percent_1)*t%span
        end if

        si%dl(:) = si%P2(:) - si%P1(:)
        si%PC(:) = si%P1(:) + si%dl(:)*(percent_c - percent_1)/(percent_2 - percent_1)
        si%zeta(:) = si%chord_c*si%dl(:)/si%ds
        si%Rroot(:) = si%PC(:) - t%root(:)
        si%af1 => t%af1
        si%af2 => t%af2
        
        !adjust start point
        if(t%side.eq.'left') then !must handle differently for opposite motion
            start(:) = si%P1(:)
        else
            start(:) = si%P2(:)
        end if

        span = span + sqrt( (si%P2(2)-si%P1(2))**2 + (si%P2(3)-si%P1(3))**2)
        t%area = t%area + si%ds
        
    end do
    t%tip(:) = start(:)
!    write(*,*) 'span: ',span
!    write(*,*) 't%area  : ',t%area
end subroutine prop_setup

!-----------------------------------------------------------------------------------------------------------
subroutine prop_write_attributes(t)
    type(prop_t) :: t

    write(*,*) '-------------------------------------------------------------'
    write(*,*) '                Propeller name : ',t%name
    write(*,*) '              Number of blades : ',t%nblades
    write(*,*) '                  Diameter (m) : [',t%diameter,']'
    write(*,*) '                     Pitch (m) : [',t%pitch,']'
    write(*,*) '                Hub Radius (m) : [',t%hubradius,']'
    write(*,*) '                Root Chord (m) : [',t%chord_1,']'
    write(*,*) '                 Tip Chord (m) : [',t%chord_2,']'
    write(*,*) '                      Rotation : [',t%rotation,']'
    write(*,*) '                  Root Airfoil : [',trim(t%af1_text),']'
    write(*,*) '                   Tip Airfoil : [',trim(t%af2_text),']'
    write(*,*) '                   Grid Points : [',t%nSec,']'
    write(*,*)

end subroutine prop_write_attributes

end module prop_m
