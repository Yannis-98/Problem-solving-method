module parametre_mod
	real,    parameter :: pi=ACOS(-1.), eps=1e-20
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

        real function fi(xi) !fonction f(x)!
	        real, intent(in) :: xi       
                fi=-cos(xi)
        end function fi

end module fonction_pb

module initialisation_mod
contains
        subroutine initialisation(nx,h,b1,b2,b3,w) !initialisation des parametre de notre equation!
        use parametre_mod
        implicit none
                integer, intent(in)  :: nx
                real,    intent(out) :: h,b1,b2,b3,w 
           	h=pi*1.0/nx
		b1=-1-h/2
		b2=2+2*h*h
		b3=-1+h/2
                w=2./(1+sin(pi*h))
        end subroutine initialisation

        subroutine initialisation_maillage(h,x) !initialisation du maillage decale!
        implicit none
                real, intent(in) :: h
                real, dimension(0:), intent(out) :: x
                integer :: nx,i
 		nx=size(x)-2
		! initialisation de x !
                do i=0, nx+1
                        x(i)=i*h-h/2
                end do		
        end subroutine initialisation_maillage

	subroutine initialisation_vecteurs(Aw, Ap, Ae, b, h, b1, b2,b3, x) !initialisation des vecteurs de notre probleme!
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



module entrees_sorties_mod !ecriture dans un fichier!
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
                        do i=1,nx
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

module ordre_methode_mod !calcul des erreurs de convergence et de consistance!
contains
	subroutine ordre_methode(nx,u, Ec, E_u, h, x)
		use fonction_pb
		implicit none
		integer :: i
		integer, intent(in) :: nx
		real, intent(in) :: h
		real, dimension(0:), intent(in) :: u, x
		real, dimension(:), intent(out) :: Ec, E_u	
		do i=1, (nx)
			Ec(i)=phi(x(i)) - u(i)
			E_u(i)=fi(x(i))+(1./(h*h))*(u(i-1)-2*u(i)+u(i+1))+(1./(2*h))*(-u(i+1)+u(i-1))-2*u(i)
		end do
		open (10, file='consistance_convergence_sor_dec.out')
		do i=1, nx
			write(10,fmt='(4e15.6)') E_u(i), Ec(i)
		end do
		close (10)

	end subroutine ordre_methode
end module ordre_methode_mod

module sor_mod
contains
	subroutine sor(s2, nx, r, b, Aw, Us, Ap, Ae, s, dU, Us_k1, w)
		use parametre_mod
		use initialisation_mod
		implicit none
		integer :: i, j
		integer, intent(in) :: nx
		real, intent(inout) :: s2, s
		real, intent(in) :: w
		real, dimension(:), intent(in) :: Ap, Aw, Ae, b
		real, dimension(:), intent(out) :: r
		real, dimension(0:), intent(out) :: Us_k1
		real, dimension(0:), intent(inout) :: Us, dU
		s=0
		s2=1e+20 !permet d'entre dans la boucle une premiere fois!
		!on initialise nos vecteurs!
		do i=0, nx+1
			Us(i)=0
		        dU(i)=0
			Us_k1(i)=0
		end do
		do i=1, nx
			r(i)=0
		end do
		do while ((s2>eps))
			Us(nx+1)=2-us(nx)
			Us(0)=-us(1)
			do j=1, nx,  1
				r(j)=b(j)-Aw(j)*Us(j-1)-Ap(j)*Us(j)-Ae(j)*Us(j+1)
				s=(r(j)*r(j))+s
			end do
			s2=sqrt(s)
			if (s2<=eps) then
			else
			do i=1,nx
				dU(i)=(w*1.)/(Ap(i))*(r(i)-Aw(i)*dU(i-1))
				Us_k1(i)=Us(i)+dU(i)
			end do
			do i=0, nx+1, 1
				Us(i)=Us_k1(i)
			end do
			end if
			s=0
		end do

	end subroutine sor
end module sor_mod

PROGRAM convergence
        use parametre_mod
        use fonction_pb
        use initialisation_mod
        use consistance_mod
        use entrees_sorties_mod
	use ordre_methode_mod
	use sor_mod
        implicit none
        integer         :: nx=1, i, j, maxi=3 !initialise le nombre de boucle!
	real, dimension(:), allocatable 	:: b, x, y, Aw, Ae, Ap
	real :: s=0,s2=1e+20, w
	REAL :: h, ti, tf
	REAL :: b1,b2,b3,NINF, const, ninf_conv, const2
	REAL, DIMENSION(:), ALLOCATABLE 	:: Lp, Lw, Ue, Ec,E_u, Us, Us_k1, r, dU 
	!Us_k1 est u^(k+1)!, et Us notre solution, dU un increment!

       do j=1,maxi
		nx=nx*10
		print *, ''
		print *, 'on affiche nx ', nx
		!on initialise!
		call initialisation(nx,h,b1,b2,b3,w)

		allocate (b(nx), x(0:nx+1), y(nx), Aw(nx), Ae(nx), Ap(nx))
		allocate (Lp(nx), Lw(nx), Ue(nx), Us(0:nx+1), dU(0:nx+1), r(nx), Us_k1(0:nx+1))
		allocate (Ec(nx),E_u(nx))


		! On initialise x pour le maillage decale !
		call initialisation_maillage(h,x)
		
		! On initialise nos vecteurs pour resoudre notre probleme !
		Call initialisation_vecteurs(Aw, Ap, Ae, b, h, b1, b2, b3, x)

		!Temps execution de SOR!
		call cpu_time(ti)
		! methode iterative, sor !
		Call sor(s2, nx, r, b, Aw, Us, Ap, Ae, s, dU, Us_k1, w)

		call cpu_time(tf)
		print *, "Temps d'execution en seconde de la methode SOR : ", tf-ti

		! ecriture de la solution Usor dans un fichier !
		if (j==maxi) then
			Call ecriture2('Solutions_SOR_dec.out',x,us, 20, '(4e15.6)')
		end if

		!Calcule des erreurs de convergence et consistance !
		Call ordre_methode(nx,Us, Ec, E_u, h, x)

		!interpretation des erreurs de convergence!
		ninf_conv=MAXVAL(ABS(Ec))
		const2=ninf_conv*nx*nx
		print *, 'La norme infini de la convergence vaut :', ninf_conv
		print *, 'La constante de la convergence vaut :', const2
		!on ecrit a la suite nos resultats de l'interpretation!
		open (unit=20, file='convergence_dec_sor.out',position='append')
		        write(unit=20,fmt='(I8, 3e15.6)') nx, h, ninf_conv, const2
		close(unit=20)
		print *, ''
		
		!interpretation des erreurs de consistance!
		NINF=MAXVAL(ABS(E_u))
		const=NINF*nx*nx
		print *, 'La norme infini de la consistance vaut :', ninf
	   	print *, 'La constante du schema avec u_LU vaut : ', const
		open (unit=20, file='consistance_dec_sor.out',position='append')
		        write(unit=20,fmt='(I8, 3e15.6)') nx, h, ninf, const
		close(unit=20)
		print *, ''
		print *, '----------------------------'

		deallocate (b, x, y, Aw, Ae, Ap, r)
		deallocate (Lp, Lw, Ue, Us, Us_k1,dU, Ec, E_u)
	end do
		!on trace la solution avec SOR!
		call execute_command_line('gnuplot -p trace_conv_dec_SOR.gnu')
	 
END PROGRAM convergence


	


