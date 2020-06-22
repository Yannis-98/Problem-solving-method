module parametre_mod
	real,    parameter :: pi=ACOS(-1.)
end module parametre_mod

module fonction_pb
contains        

        real function phi(xi) !calcul de la solution exacte !
        use parametre_mod 
        implicit none
                real, intent(in) :: xi        
                real ::  alpha
                alpha=(7-3*exp(2*pi))/(10.0*(exp(-pi)-exp(2*pi)))    
                phi=alpha*(exp(-xi)-exp(2*xi))+(3.0/10)*exp(2*xi)-(3/10.0)*cos(xi)-(1/10.0)*sin(xi)
        end function phi

        real function fi(xi) !calcul de f(x)!
	        real, intent(in) :: xi       
                fi=-cos(xi)
        end function fi

end module fonction_pb

module initialisation_mod
contains

        subroutine initialisation(nx,h,a1,a2,a3) !initialisation des parametres de notre probleme !
        use parametre_mod
        implicit none
                integer, intent(in)  :: nx
                real,    intent(out) :: h,a1,a2,a3
                h=pi*1.0/nx
                a1=1+h/2.0 
                a2=2+2*h*h
                a3=-1+h/2.0
        end subroutine initialisation

        subroutine initialisation_maillage(h,x)
        implicit none
                real, intent(in) :: h
                real, dimension(0:), intent(out) :: x
                integer :: nx,i
                nx=size(x)-1
	! initialisation de x !
                do i=0, nx
                        x(i)=i*h
                end do
        end subroutine initialisation_maillage

	subroutine initialisation_vecteurs(Aw, Ap, Ae, b, h, a1, a2,a3, x)
		use fonction_pb
		implicit none
		integer :: i, nx
		real, dimension(:), allocatable 	:: x, b,Aw, Ae, Ap
		REAL :: a1,a2,a3, h
		!Initialisation de nos vecteur du probleme !
		nx=size(x)-1
		DO i=1, (nx-1), 1
			Aw(i)=-a1
			Ap(i)=a2
			Ae(i)=a3
			b(i)=h*h*(fi(x(i)))
		END DO
		b(nx-1)=b(nx-1)-a3
		Aw(1)=0
		Ae(nx-1)=0
	end subroutine initialisation_vecteurs

end module initialisation_mod


module entrees_sorties_mod
contains

        subroutine ecriture(nom_fichier,v1,v2, n, format1) !Ecriture d'un fichier avec le nom en parametre de la fonction !
        use fonction_pb
        implicit none
                character(len=*), intent(in) :: nom_fichier, format1
                real, dimension(0:), intent(in) :: v1
                real, dimension(:), intent(in) :: v2
		integer :: i,nx, n
                nx=size(v1)-1
                open(unit=n,file=nom_fichier) 
                        do i=1,nx-1
                                write(unit=n,fmt=format1) v1(i),fi(v1(i)),phi(v1(i)),v2(i)
                        enddo
                close(unit=n)
        end subroutine ecriture

	subroutine ecriture2(nom_fichier,v1,v2, n, format1) !Ecriture d'un fichier avec le nom en parametre de la fonction !
        use fonction_pb
        implicit none
                character(len=*), intent(in) :: nom_fichier, format1
                real, dimension(0:), intent(in) :: v1
                real, dimension(0:), intent(in) :: v2
		integer :: i,nx, n
                nx=size(v1)-1
                open(unit=n,file=nom_fichier)
                        do i=0,nx
                                write(unit=n,fmt=format1) v1(i),v2(i)
                        enddo
                close(unit=n)
        end subroutine ecriture2

end module entrees_sorties_mod

module LU
contains

	subroutine initialisation_LU(Lp, Lw, Ue, Ap, Ae, Aw, nx) !initialisation des vecteurs de notre matrice LU !
		implicit none
		integer, intent(in) :: nx
		integer :: i
		real, dimension(:), intent(in) 	::  Aw, Ae, Ap
		real, dimension(:), intent(out) :: Lp, Lw, Ue
		Lp(1)=Ap(1)
		Ue(1)=Ae(1)*1.0/Lp(1)
		DO i=2, (nx-2), 1
			Lw(i)=Aw(i)
			Lp(i)=Ap(i)-Lw(i)*Ue(i-1)
			Ue(i)=Ae(i)*1.0/Lp(i)
		end do
		Lw(nx-1)=Aw(nx-1)
		Lp(nx-1)=Ap(nx-1)-Lw(nx-1)*Ue(nx-2)
	end subroutine initialisation_LU

	subroutine descente_LU(y,b,Lp,Lw, nx) !descente LU !
		implicit none
		integer, intent(in) :: nx
		integer :: i
		real, dimension(:), intent(in) 	:: b, Lp, Lw
		real, dimension(:), intent(out) :: y
		y(1)=b(1)*1.0/Lp(1)
		DO i=2, (nx-1), 1
			y(i)=(b(i)-Lw(i)*y(i-1))/Lp(i)
		end do
	end subroutine descente_LU

	subroutine remontee_LU(u,y,Ue, nx) !remontee LU !
		implicit none
		integer, intent(in) :: nx
		integer :: i
		real, dimension(:), intent(in) 	::  Ue,y
		real, dimension(0:), intent(out) :: u
		u(nx-1)=y(nx-1)
		DO i=(nx-2), 1, -1
			u(i)=y(i)-Ue(i)*u(i+1)
		end do
		u(0)=0
		u(nx)=1
	end subroutine remontee_LU

end module LU

module ordre_methode_mod
contains

	subroutine ordre_methode(nx,u, Ec, E_u, h, x)
		use fonction_pb
		implicit none
		integer :: i
		integer, intent(in) :: nx
		real, intent(in) :: h
		real, dimension(0:), intent(in) :: u, x
		real, dimension(:), intent(out) :: Ec, E_u	
		do i=1, (nx-1) !calcul des erreurs de convergence et consistance!
			Ec(i)=phi(x(i)) - u(i)
			E_u(i)=fi(x(i))+(1./(h*h))*(u(i-1)-2*u(i)+u(i+1))+(1./(2*h))*(-u(i+1)+u(i-1))-2*u(i)
		end do
		open (10, file='consistance_convergence_dec_LU.out')  !stockage dans un fichier!
		do i=1, nx-1
			write(10,fmt='(4e15.6)') E_u(i), Ec(i)
		end do
		close (10)
	end subroutine ordre_methode

end module ordre_methode_mod


PROGRAM convergence
        use parametre_mod
        use fonction_pb
        use initialisation_mod
        use entrees_sorties_mod
	use LU
	use ordre_methode_mod

        implicit none
 
        integer         :: nx=1, i, j, maxi=6 !nombre d'iteration!
	real, dimension(:), allocatable 	:: b, x, y, Aw, Ae, Ap
	REAL :: h, ti, tf
	REAL :: a1,a2,a3,NINF, const, ninf_conv, const2

	REAL, DIMENSION(:), ALLOCATABLE 	:: Lp, Lw, Ue, u, Ec,E_u
	!Ec est erreur de convergence et E_u erreur de constistance avec u !
	!Ue est pour LU, u est la solution !

        do j=1,maxi
		nx=nx*10
		print *, 'Pour nx=', nx

		call initialisation(nx,h,a1,a2,a3)

		allocate (b(nx-1), x(0:nx), y(nx-1), Aw(nx-1), Ae(nx-1), Ap(nx-1))
		allocate (Lp(nx-1), Lw(nx-1), Ue(nx-1), u(0:nx))
		allocate (Ec(nx-1),E_u(nx-1))
		call initialisation_maillage(h,x)
		Call initialisation_vecteurs(Aw, Ap, Ae, b, h, a1, a2, a3, x)
	!Résolution LU !
		call cpu_time(ti) !temps d'execution de LU!

		!Initialisation LU !
		Call initialisation_LU(Lp, Lw, Ue, Ap, Ae, Aw, nx)
		! Déscente !
		Call descente_LU(y,b,Lp,Lw, nx)
		! Remontée !
		Call remontee_LU(u,y,Ue, nx)

		! ecriture de la solution u dans un fichier pour hmax!
		if (j==maxi) then 
			Call ecriture2('Solutions_LU.out',x,u, 20, '(4e15.6)')
			print *, 'Fichier consistance_cla_LU.out et convergence_cla_LU.out et Solutions_LU.out OK '
		end if

		call cpu_time(tf) !temps execution lU!
		print *, "Temps d'execution en seconde de la methode LU : ", tf-ti
		print *, ''

		!Calcul des Erreurs de convergence et consistance!
		Call ordre_methode(nx,u, Ec, E_u, h, x)

		!Erreurs de convergence!
		ninf_conv=MAXVAL(ABS(Ec))
		const2=ninf_conv*nx*nx
		print *, 'La norme infini de la convergence vaut :', ninf_conv
		print *, 'La constante de la convergence vaut :', const2
		print *, ''
		open (unit=20, file='convergence_cla_LU.out',position='append')
		        write(unit=20,fmt='(I8, 3e15.6)') nx, h, ninf_conv, const2
		close(unit=20)

		!Erreurs de consistance!
		NINF=MAXVAL(ABS(E_u))
		const=NINF*nx*nx
		print *, 'La norme infini de la consistance vaut :', ninf
	    	print *, 'La constante du schema avec u_LU vaut : ', const
		open (unit=20, file='consistance_cla_LU.out',position='append')
		        write(unit=20,fmt='(I8, 3e15.6)') nx, h, ninf, const
		close(unit=20)
		print *, ''

		deallocate (b, x, y, Aw, Ae, Ap)
		deallocate (Lp, Lw, Ue, u)
		deallocate (Ec,E_u)
		print *, '------------------------'
	end do
	call execute_command_line('gnuplot -p trace_conv_clas_LU.gnu') !trace la solution avec le dernier h calculé!
 
END PROGRAM convergence


	


