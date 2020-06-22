module parametre_mod
	real,    parameter :: pi=ACOS(-1.)
end module parametre_mod

module initialisation_mod
contains
        subroutine initialisation(nx,h,a1,a2,a3,w)
        use parametre_mod
        implicit none
                integer, intent(in)  :: nx
                real,    intent(out) :: h,a1,a2,a3,w 
! initialisation de la matrice du probleme
                h=pi*1.0/nx
                a1=1+h/2.0 
                a2=2+2*h*h
                a3=-1+h/2.0
                w=2./(1+sin(pi*h))
        end subroutine initialisation

        subroutine initialisation_maillage(h,x)
        implicit none
                real, intent(in) :: h
                real, dimension(0:), intent(out) :: x
                integer :: nx,i
                nx=size(x)-1
	! initialisation de x pour le maillage classique !
                do i=0, nx
                        x(i)=i*h
                end do

        end subroutine initialisation_maillage
end module initialisation_mod

module fonction_pb
contains        
        real function phi(xi) ! Fonction permettant de calculer la solution theorique !
        use parametre_mod 
        implicit none
                real, intent(in) :: xi
                real ::  alpha
                alpha=(7-3*exp(2*pi))/(10.0*(exp(-pi)-exp(2*pi)))
                phi=alpha*(exp(-xi)-exp(2*xi))+(3.0/10)*exp(2*xi)-(3/10.0)*cos(xi)-(1/10.0)*sin(xi)

        end function phi
        real function fi(xi) ! Fonction permettant de calculer la fonction f(x) !
	        real, intent(in) :: xi
                fi=-cos(xi)
        end function fi

end module fonction_pb

module entrees_sorties_mod
contains
        subroutine ecriture(nom_fichier,v1,v2) !Permet d'ecrire sur un fichier !
        use fonction_pb
        implicit none
                character(len=*), intent(in) :: nom_fichier
                real, dimension(0:), intent(in) :: v1
                real, dimension(:), intent(in) :: v2
                integer :: i,nx
                nx=size(v1)-1
                open(unit=10,file=nom_fichier) ! nom du fichier !
                        do i=1,nx-1
                                write(unit=10,fmt='(4e15.6)') v1(i),fi(v1(i)),phi(v1(i)),v2(i)
                        enddo
                close(unit=10)
        end subroutine ecriture
end module entrees_sorties_mod

module consistance_mod
contains
        subroutine erreur(h,x,E) ! Calcul de la consistance !
        use fonction_pb
        implicit none
                real, intent(in) :: h
                real, dimension(0:), intent(in) :: x
                real, dimension(:), intent(out) :: E !erreur de consistance !
                integer :: nx, i
                nx = size(x)-1
                do i=1, (nx-1), 1
                        E(i)=fi(x(i))+ &
                           (1./(h*h))*(phi(x(i-1))-2*phi(x(i))+phi(x(i+1)))+(1./(2*h))*(-phi(x(i+1))+phi(x(i-1)))-2*phi(x(i))
                end do
        end subroutine erreur
end module consistance_mod

PROGRAM consistance
        use parametre_mod
        use fonction_pb
        use initialisation_mod
        use consistance_mod
        use entrees_sorties_mod
        implicit none
        integer         :: nx=1, i, k, j
	real, dimension(:), allocatable 	:: x, E
	REAL :: h
	REAL :: a1,a2,a3,NINF, const,w
       	do j=1,6
		nx=nx*10
		allocate (x(0:nx))
		allocate (E(nx-1))
		call initialisation(nx,h,a1,a2,a3,w)
		call initialisation_maillage(h,x)
		call erreur(h,x,E) ! calcul de l'erreur de consistance !
		call ecriture("consistance.out",x,E)

		NINF=MAXVAL(ABS(E)) ! calcul de la consistance !
		const=NINF*nx*nx

		print *, ' '
		print *, 'pour nx =', nx
		print *, ' '
		print *, 'La constante du schema vaut : ', const
		open (unit=20, file='ordre.out',position='append') ! on ecrit a la suite du fichier la constante de consistance
		        write(unit=20,fmt='(I8, 3e15.6)') nx, h, ninf, const
		close(unit=20)
		print *, ' '
		print *, '-------------------'
		print *, ' '
		deallocate(x,E)

	end do
		call execute_command_line('gnuplot -p trace_ordre_h_ninf.gnu') !tracer de ninf en fonction de h via script
		call execute_command_line('gnuplot -p trace_ordre_h_const.gnu') !tracer de const en fonction de h via script
	
!------------------------------!
END PROGRAM consistance

