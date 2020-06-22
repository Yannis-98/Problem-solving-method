module parametre_mod
	real,    parameter :: pi=ACOS(-1.)
end module parametre_mod

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

        real function fi(xi) !fonction f(x) !
	        real, intent(in) :: xi
                fi=-cos(xi)
        end function fi

end module fonction_pb

module initialisation_mod
contains
        subroutine initialisation(nx,h,b1,b2,b3) !parametrage des equations!
        use parametre_mod
        implicit none
                integer, intent(in)  :: nx
                real,    intent(out) :: h,b1,b2,b3 
           	h=pi*1.0/nx
		b1=-1-h/2
		b2=2+2*h*h
		b3=-1+h/2
        end subroutine initialisation

        subroutine initialisation_maillage(h,x) !initialisation maillage decale!
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

	subroutine initialisation_vecteurs(Aw, Ap, Ae, b, h, b1, b2,b3, x) !Initialisation de nos vecteur via dvp analytique!
		use fonction_pb
		implicit none
		integer :: i, nx
		real, dimension(:), allocatable 	:: x, b,Aw, Ae, Ap
		REAL :: b1,b2,b3, h
		nx=size(x)-2
		DO i=2, (nx-1), 1
			Aw(i)=b1
			Ap(i)=b2
			Ae(i)=b3
			b(i)=h*h*(fi(x(i)))
		END DO
		b(1)=h*h*(fi(x(1)))
		b(nx)=h*h*(fi(x(nx)))-2*b3
		Ap(1)=b2-b1
		Ap(nx)=b2-b3
		Aw(1)=0
		Aw(nx)=b1
		Ae(1)=b3
		Ae(nx)=0
	end subroutine initialisation_vecteurs

end module initialisation_mod

module entrees_sorties_mod !ecriture des fichiers!
contains
        subroutine ecriture(nom_fichier,v1,v2, n, format1)
        use fonction_pb
        implicit none
                character(len=*), intent(in) :: nom_fichier, format1
                real, dimension(0:), intent(in) :: v1
                real, dimension(:), intent(in) :: v2
		integer :: i,nx, n
                nx=size(v1)-2
                open(unit=n,file=nom_fichier)
                        do i=1,nx-1
                                write(unit=n,fmt=format1) v1(i),fi(v1(i)),phi(v1(i)),v2(i)
                        enddo
                close(unit=n)
        end subroutine ecriture

	subroutine ecriture2(nom_fichier,v1,v2, n, format1)
        use fonction_pb
        implicit none
                character(len=*), intent(in) :: nom_fichier, format1
                real, dimension(0:), intent(in) :: v1
                real, dimension(0:), intent(in) :: v2
		integer :: i,nx, n
                nx=size(v1)-2
                open(unit=n,file=nom_fichier)
                        do i=0,nx+1
                                write(unit=n,fmt=format1) v1(i),v2(i)
                        enddo
                close(unit=n)
        end subroutine ecriture2

end module entrees_sorties_mod

module LU
contains
	subroutine initialisation_LU(Lp, Lw, Ue, Ap, Ae, Aw, nx) !initialisation des matrice L et U avec des vecteurs!
		implicit none
		integer, intent(in) :: nx
		integer :: i
		real, dimension(:), intent(in) 	::  Aw, Ae, Ap
		real, dimension(:), intent(out) :: Lp, Lw, Ue
		Lp(1)=Ap(1)
		Ue(1)=Ae(1)*1.0/Lp(1)
		DO i=2, (nx-1), 1
			Lw(i)=Aw(i)
			Lp(i)=Ap(i)-Lw(i)*Ue(i-1)
			Ue(i)=Ae(i)*1.0/Lp(i)
		end do
		Lw(nx)=Aw(nx)
		Lp(nx)=Ap(nx)-Lw(nx)*Ue(nx-1)
	end subroutine initialisation_LU
	subroutine descente_LU(y,b,Lp,Lw, nx) !descente LU !
		implicit none
		integer, intent(in) :: nx
		integer :: i
		real, dimension(:), intent(in) 	:: b, Lp, Lw
		real, dimension(:), intent(out) :: y
		y(1)=b(1)*1.0/Lp(1)
		DO i=2, (nx), 1
			y(i)=(b(i)-Lw(i)*y(i-1))/Lp(i)
		end do
	end subroutine descente_LU
	subroutine remontee_LU(u,y,Ue, nx) !remontee lU !
		implicit none
		integer, intent(in) :: nx
		integer :: i
		real, dimension(:), intent(in) 	::  Ue,y
		real, dimension(0:), intent(out) :: u
		u(nx)=y(nx)
		DO i=(nx-1), 1, -1
			u(i)=y(i)-Ue(i)*u(i+1)
		end do
		!extrapolation pour maillage decale!
		u(0)=-u(1)
		u(nx+1)=2-u(nx)
	end subroutine remontee_LU

end module LU

module ordre_methode_mod
contains
	subroutine ordre_methode(nx,u, Ec, E_u, h, x) !calcul erreurs de convergence et consistance!
		use fonction_pb
		implicit none
		integer :: i
		integer, intent(in) :: nx
		real, intent(in) :: h
		real, dimension(0:), intent(in) :: u, x
		real, dimension(:), intent(out) :: Ec, E_u	
		do i=1, nx
			Ec(i)=phi(x(i)) - u(i)
			E_u(i)=fi(x(i))+(1./(h*h))*(u(i-1)-2*u(i)+u(i+1))+(1./(2*h))*(-u(i+1)+u(i-1))-2*u(i)
		end do
		open (10, file='consistance_convergence_dec.out')
		do i=1, nx
			write(10,fmt='(4e15.6)') E_u(i), Ec(i)
		end do
		close (10)
	end subroutine ordre_methode
end module ordre_methode_mod



PROGRAM convergence
        use parametre_mod
        use fonction_pb
        use initialisation_mod
        use consistance_mod
        use entrees_sorties_mod
	use ordre_methode_mod
	use LU
        implicit none
 
        integer         :: nx=1, i, j, maxi=6 !nb iteration !
	real, dimension(:), allocatable 	:: b, x, y, Aw, Ae, Ap
	REAL :: h, ti, tf
	REAL :: b1,b2,b3,NINF, const, ninf_conv, const2
	REAL, DIMENSION(:), ALLOCATABLE 	:: Lp, Lw, Ue, u, Ec,E_u

        do j=1,maxi
		nx=nx*10
		print *, 'Pour nx =', nx
		print *, ' '

		call initialisation(nx,h,b1,b2,b3)

		allocate (b(nx), x(0:nx+1), y(nx), Aw(nx), Ae(nx), Ap(nx))
		allocate (Lp(nx), Lw(nx), Ue(nx), u(0:nx+1))
		allocate (Ec(nx),E_u(nx))
		!initialisation!
		call initialisation_maillage(h,x)
		Call initialisation_vecteurs(Aw, Ap, Ae, b, h, b1, b2, b3, x)

		! On mesure le temps d execution de LU!
		call cpu_time(ti)

		!Initialisation LU !
		Call initialisation_LU(Lp, Lw, Ue, Ap, Ae, Aw, nx)
		!Résolution LU !
		! Déscente !
		Call descente_LU(y,b,Lp,Lw, nx)
		! Remontée !
		Call remontee_LU(u,y,Ue, nx)

		call cpu_time(tf)
		print *, "Temps d'execution en seconde de la methode LU : ", tf-ti

		! ecriture de u dans un fichier pour la plus grande valeur de h!
		if (j==maxi) then
			Call ecriture2('Solutions_LU_dec.out',x,u, 20, '(4e15.6)')
			call execute_command_line('gnuplot -p trace_conv_dec_LU.gnu')
		end if

		!Calcul des erreurs de convergence et consistance !
		Call ordre_methode(nx,u, Ec, E_u, h, x)
		print *, ''
		
		!interpretation des erreurs de convergence !
		ninf_conv=MAXVAL(ABS(Ec))
		const2=ninf_conv*nx*nx
		print *, 'La norme infini de la convergence vaut :', ninf_conv
		print *, 'La constante de la convergence vaut :', const2
		print *, ''
		open (unit=20, file='convergence_dec_LU.out',position='append')
		        write(unit=20,fmt='(I8, 3e15.6)') nx, h, ninf_conv, const2
		close(unit=20)

		!interpretation des erreurs de consistance !
		NINF=MAXVAL(ABS(E_u))
		const=NINF*nx*nx
		print *, 'La norme infini de la consistance vaut :', ninf
	   	print *, 'La constante du schema avec u_LU vaut : ', const
		open (unit=20, file='consistance_dec_LU.out',position='append')
		        write(unit=20,fmt='(I8, 3e15.6)') nx, h, ninf, const
		close(unit=20)
		print *, ''
		print *, '---------------'
		print *, ''

		deallocate (b, x, y, Aw, Ae, Ap)
		deallocate (Lp, Lw, Ue, u)
		deallocate (Ec,E_u)

	end do

 
END PROGRAM convergence


	


