module parametre_mod
	real,    parameter :: pi=ACOS(-1.)
end module parametre_mod

module initialisation_mod
contains
        subroutine initialisation(nx,h,b1,b2,b3,w) !initialisation des parametre de l'equation discrete!
        use parametre_mod
        implicit none
                integer, intent(in)  :: nx
                real,    intent(out) :: h,b1,b2,b3,w 
           	h=pi*1.0/nx
		b1=-1
		b2=2+2*h*h-h
		b3=h-1
                w=2./(1+sin(pi*h))
        end subroutine initialisation

        subroutine initialisation_maillage(h,x) !initialisation du maillage decale!
	use parametre_mod
        implicit none
                real, intent(in) :: h
                real, dimension(0:), intent(out) :: x
                integer :: nx,i
                nx=size(x)-2
                do i=0, nx+1
                        x(i)=i*h-h/2
                end do
        end subroutine initialisation_maillage

end module initialisation_mod

module fonction_pb
contains        
        real function phi(xi) !solution exacte!
        use parametre_mod 
        implicit none
                real, intent(in) :: xi
                real ::  alpha        
                alpha=(7-3*exp(2*pi))/(10.0*(exp(-pi)-exp(2*pi)))        
                phi=alpha*(exp(-xi)-exp(2*xi))+(3.0/10)*exp(2*xi)-(3/10.0)*cos(xi)-(1/10.0)*sin(xi)
        end function phi

        real function fi(xi) !fonction f(x)!
	        real, intent(in) :: xi      
                fi=-cos(xi)
        end function fi

end module fonction_pb

module entrees_sorties_mod !ecriture des fichiers!
contains
        subroutine ecriture(nom_fichier,v1,v2) 
        use fonction_pb
        implicit none
                character(len=*), intent(in) :: nom_fichier
                real, dimension(0:), intent(in) :: v1
                real, dimension(:), intent(in) :: v2
                integer :: i,nx
                nx=size(v1)-2
                open(unit=40,file=nom_fichier)
                        do i=1,nx
                                write(unit=40,fmt='(4e15.6)') v1(i),fi(v1(i)),phi(v1(i)),v2(i)
                        enddo
                close(unit=40)
        end subroutine ecriture
end module entrees_sorties_mod

module consistance_mod
contains
        subroutine erreur(h,x,E) !calcul de l'erreur de consistance du schema decale!
        use fonction_pb
        implicit none
                real, intent(in) :: h
                real, dimension(0:), intent(in) :: x
                real, dimension(:), intent(out) :: E 
                integer :: nx, i
                nx = size(x)-2
                do i=1, nx, 1
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
        integer         :: nx=1, i, j
	real, dimension(:), allocatable 	:: x, E
	REAL :: h
	REAL :: b1,b2,b3,NINF, const,w 
        do j=1,6
		nx=nx*10 !evolution de de nx!
		print *, 'Pour nx =', nx

		call initialisation(nx,h,b1,b2,b3,w) 

		allocate (x(0:nx+1))
		allocate (E(nx))

		call initialisation_maillage(h,x)
		!calcul de l'erreur de consistance et sauvegarde !
		call erreur(h,x,E)
		call ecriture("difference.out",x,E)

		!interpretation de l'erreur de consistance!
		NINF=MAXVAL(ABS(E))
		const=NINF*(nx)*(nx)
		print *, ' '
		print *, 'La constante du schema decale vaut : ', const
		open (unit=20, file='ordre_dec.out',position='append') !stockage a la suite !
		        write(unit=20,fmt='(I8, 3e15.6)') nx, h, ninf, const
		close(unit=20)
	

		deallocate (x)
		deallocate (E)
		print *, '--------------------'
		print *, ' '
	end do

	call execute_command_line('gnuplot -p trace_ordre2_h_ninf.gnu') !on trace ninf en fonction de h pour la consistance!
	call execute_command_line('gnuplot -p trace_ordre2_h_const.gnu') !on trace const en fonction de h pour consistance !

END PROGRAM consistance

