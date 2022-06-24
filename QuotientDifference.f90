Program QuotientDifference

!USE METODOSQD
!1- Ingresa el grado del polinomio
!2- Ingresa los coeficientes de forma ordenada
!3- Comienza QD
    !i-Ejecuta primeras iteraciones

IMPLICIT NONE
!----------------------
integer ,parameter:: grado=5
!---------------------------
integer i ,Metodo,maxiter
real(8), dimension(0:grado) :: VecCoeficientes !Hay un coef  mas que el grado
complex(8), dimension(1:grado) :: VecRaices 
complex(8), dimension(1:grado) :: VecRaicesCOMPLEJAS
real(8) tol,diferencia

! COMIENZA EL PROGRAMA
!-----------

!-------------
CALL LeeCoeficientes(VecCoeficientes,Grado)
!----------ERROR--------!
!write(*,*)'Ingrese la cota maxima de error'
!read(*,*)tol
tol=0.00001
!------Iteraciones-----!
!write(*,*)"Ingrese maxima cantidad de iteraciones" 
!read(*,*)maxiter
maxiter=10


Metodo=0

!CALL PreProcesamiento(VecCoeficientes,grado,Metodo,diferencia)


do while (CoeficienteCero(VecCoeficientes,Grado+1) .AND. Metodo<4) 

	Metodo=Metodo+1
	!CALL PreProcesamiento(VecCoeficientes,grado,Metodo,Diferencia) !hay 3 metodos para corregir 0s, haciendo diferentes reempazos s=x+d s+cx s=1/x
	                                                                  !ademas puede bajar de grado si el ultimo termino es 0
	write(*,*)"Entro y no deberia"
enddo

if (CoeficienteCero(VecCoeficientes,Grado+1)) then !si el metodo es 4 enotnces no fue posible alterar el polinomio

  write(*,*)"Imposible obtener raices del polinomio : No fue posible modificarlo para que sus coeficientes fueran distintos de 0"
  
 else
	!select case (grado)
	!case (1) !grado=0
	!	write(*,*)"Error grado=0"
	!case (2)
	!	write(*,*)"Raiz: ",VecCoeficientes(1)/VecCoeficientes(2)
	!case (3)
		!call Resolvente(VecCoeficientes)!Falta desarrollar
	!case default
		call QuotientDifferenceSub(VecCoeficientes,VecRaices,grado,tol,maxiter)
		!call CorrigeSolucion(VecRaices,grado,Metodo,diferencia)
		!call MuestraRaices(VecRaices,VecRaicesCOMPLEJAS,grado)
	!end select
endif




CONTAINS
!----------------------------------!
subroutine MuestraRaices(VecRaices,grado)
complex(8), intent(in) :: VecRaices(:)
integer i,grado

do i=1, grado
		write(*,*)"Raiz ",i,": ",VecRaices(i)
enddo



end subroutine
!----------------------------------!
subroutine PreProcesamiento(VecCoeficientes,Grado,Metodo,diferencia)
real(8), intent(inout) :: VecCoeficientes(:)
real(8), intent(out) :: diferencia
integer, intent(inout) ::Metodo,Grado

diferencia = 1 ! variable para el metodo 1 x-d  darle algun valor 
select case(Metodo)
	case (0) !Caso 0 en ultimo lugar
		do while (VecCoeficientes(grado)==0)
			grado=grado-1
		enddo
	case (1)
		call Metodo1(VecCoeficientes,diferencia,Grado)
	case (2)
	case (3)
end select

endsubroutine
!----------------------------------!
subroutine Metodo1(VecCoeficientes,d,Grado)
real(8) ,intent(inout) :: VecCoeficientes(:),d
integer, intent(in) :: Grado
integer i,j,k

do i=2,grado+1

 k=i
 if(VecCoeficientes(i)/=0) then
   do j=1,i-1
    VecCoeficientes(j) = VecCoeficientes(j) + Combinatoria(i-1,k-1)*(VecCoeficientes(i)**(i-k))*(-1*d)**(k-1)
    k=k-1
   end do
 
  VecCoeficientes(j+1) = (VecCoeficientes(i)**(i-1))
 end if

end do

endsubroutine Metodo1
!----------------------------------!
function Combinatoria(m,n)
integer, intent(in) :: m,n
real(8) Combinatoria

Combinatoria=Factorial(m)/((Factorial(n)*Factorial(m-n))*1.0) !el resultado de la operacion es un entero asi q *1.0 para q sea real

endfunction Combinatoria
!----------------------------------!
function Factorial(x)
integer, intent(in) :: x
integer Factorial

if (x==0) then
	Factorial=1
else
	Factorial=1
	do i=1,x
		Factorial=Factorial*i
	enddo
end if	

end function Factorial
!----------------------------------!

function CoeficienteCero(VecCoeficientes,Grado)
logical CoeficienteCero
real(8) VecCoeficientes(:)
integer Grado,i

do while ((VecCoeficientes(i)/=0) .AND. (i<Grado))
	i=i+1
enddo

CoeficienteCero=((Grado==i).AND.(VecCoeficientes(grado)/=0))

endfunction

!----------------------------------!
SUBROUTINE LeeCoeficientes(VecCoeficientes,Grado)
integer i,Grado
real(8) :: VecCoeficientes(:)

open(2,FILE="CoeficientesPolinomio.txt", action="read")


write(*,'(A,I2,A)',advance='NO')"Polinomio de grado ",grado," : "
do i=1,grado+1 ,1
	read(2,*)VecCoeficientes(i)
	WRITE(*,'(F7.2,A,I2,A)',ADVANCE='NO')VecCoeficientes(i),"x^",i-1," + " !'(F4.2,A,I,A)'
	!VecCoeficientes(i)=VecCoeficientes(i)*((1/5.)**i)
enddo

write(*,*)

END SUBROUTINE
!----------------------------------!

subroutine QuotientDifferenceSub(VecCoeficientes,VecRaices,grado,tol,maxiter)
!parametros 
real(8) ,intent(inout)::VecCoeficientes(:)

complex(8), dimension(1:grado), intent(inout) :: VecRaices
integer ,intent(in):: grado,maxiter
real(8) ,intent(in)::tol
!var locales
complex(8), dimension(1:grado) :: VecQ_ANT
real(8) ,dimension(1:grado) :: VecQ
integer i,iter
real(8) error
real(8),dimension(0:grado) :: VecError
real(8) CoefV,CoefU

complex(8) RaizComp1,RaizComp2
!VECError va de 0 a grado
!VEcRaiz va de 1 a grado
!VecCoeficientes va de 1 a grado

   
!COMIENZA QD
CALL CondicionesIniciales(VecCoeficientes,VecQ,VecError,grado) 

iter=1
error=2*tol 
write(*,*)"iter",iter

do while (error>tol .AND. maxiter>iter)
write(*,*)"-------------------------"
	VecQ_ANT=VecQ
    do i=1, grado   !Calculo la fila de entera en cada iteración y el grado +1 porque el vector comiena en 1 , el el coef de maximo grado esta een grado +1
       VecQ(i)=VecError(i)-VecError(i-1)+VecQ(i)     
       
    enddo
    
    
    Do i=1, grado-1      !Como e(0)=e(n)=0 siempre, uso esos límites en el ciclo
       VecError(i)=VecError(i)*(VecQ(i+1)/VecQ(i))   ! Hay un error de tipeo en el pdf que figura q(i-1) pero en la tabla utiliza q(i+1), sino no da
    end do
    

     iter=iter+1
     error=maxval(abs(VecError))
     write(*,*)"iter",iter
     write(*,*)"Raices : ",VecQ
     write(*,*)"Errores : ",VecError
end do 

		do i=0,grado-1
			if (abs(VecError(i))>tol) then
				CoefU=VecQ(i)+VecQ(i+1)
				CoefV=(-1)*VecQ_ANT(i)*VecQ(i+1) !No se si va el -1 pero lo vi asi
				write(*,*)CoefU
				write(*,*)CoefV
				CALL MetodoBairstow (VecCoeficientes,grado,tol,CoefU,CoefV)
				CALL Resolvente(-CoefU, -CoefV, RaizComp1, RaizComp2)
				VecRaices(i)=RaizComp1
				VecRaices(i+1)=RaizComp2
			else
				VecRaices(i+1)=DCMPLX(VecQ(i+1), 0.)
			endif
		enddo
	
write(*,*)VecRaices
end subroutine

!----------------------------------!

subroutine CondicionesIniciales(VecCoeficientes,VecRaices,VecError,grado) !Esta Bien
real(8), intent(inout)::VecCoeficientes(0:),VecRaices(1:)
integer ,intent(in):: grado
real(8), intent(inout) :: VecError(0:) !NO SACAR ESTE 0, SI LO SACAS REVIENTA EL PROGRAMA
integer i

VecRaices=0.
VecRaices(1)=-1.*(VecCoeficientes(grado-1))/(VecCoeficientes(grado))

do i=1, grado-1

      VecError(i)=VecCoeficientes(grado-i-1)/VecCoeficientes(grado-i)
      !write(*,*)VecCoeficientes(grado+1-i-1),VecCoeficientes(grado+1-i),"  ",grado-i-1,"i=",i," grado=",grado,VecError(i)
      
end do

VecError(0)=0.
VecError(grado)=0.

write(*,*)"Raices : ",VecRaices
     write(*,*)"Errores : ",VecError
endsubroutine

!----------------------------------!

subroutine MetodoBairstow (VecCoeficientes,grado,tol,CoefU,CoefV)
real(8), intent(in) :: VecCoeficientes(0:)
real(8), intent(inout) :: CoefU,CoefV
real(8), intent(in) :: tol
integer, intent(in) ::grado
!var locales

real(8) Q, Q_ant1, Q_ant2 ,P, P_ant1, P_ant2, P_ant3, IncrementoH, IncrementoK, error
integer i

error=2*tol
do while (error >= tol)
	Q_ant1 = 0.
    Q_ant2 = 0.
    P_ant1 = 0.
    P_ant2 = 0.
	do i=grado, 1, -1
	
		Q = VecCoeficientes(i) + CoefU*Q_ant1 + CoefV*Q_ant2
		Q_ant2 = Q_ant1
		Q_ant1 = Q
		
		P = Q + CoefU*P_ant1 + CoefV*P_ant2
		P_ant3 = P_ant2
        P_ant2 = P_ant1
        P_ant1= P
	enddo
	q = VecCoeficientes(0) + CoefU* Q_ant1 + CoefV* Q_ant2
    !write(*,*)"Q ",Q,"Q_ant1 ",Q_ant1,"Q_ant2",Q_ant2,"P: ",P,"P_ant1",P_ant1,"P_ant2",P_ant2,"P_ant3",P_ant3
	IncrementoH = (Q * P_ant3 - Q_ant1 * P_ant2) / (P_ant2**2. - P_ant1 * P_ant3)
	IncrementoK = (Q_ant1 * P_ant1 - Q * P_ant2) / (P_ant2**2. - P_ant1 * P_ant3)
	!write(*,*)"Incremento K: ",IncrementoK
	CoefU = CoefU +IncrementoH
    CoefV = CoefV +IncrementoK

	if(abs(Q) > abs(Q_ant1))then
		error = abs(Q)
	else
		error= abs(Q_ant1)
	endif
	
enddo

end subroutine

!----------------------------------!

SUBROUTINE Resolvente(b, c, c1, c2)

REAL(8) a
REAL(8), INTENT(IN) :: b, c
COMPLEX(8) Discriminante 
COMPLEX(8), INTENT(OUT) :: c1, c2

    a = 1.
    Discriminante = b**2 - 4. * a * c
    c1 = (-b + sqrt(Discriminante)) / (2. * a)
    c2 = (-b - sqrt(Discriminante)) / (2. * a)

END SUBROUTINE Resolvente



!----------------------------------!
End program
