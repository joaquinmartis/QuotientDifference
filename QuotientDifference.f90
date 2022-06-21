Program QuotientDifference

!1- Ingresa el grado del polinomio
!2- Ingresa los coeficientes de forma ordenada
!3- Comienza QD
    !i-Ejecuta primeras iteraciones

IMPLICIT NONE

integer:: Grado, i ,Metodo
real(8), allocatable, dimension(:) :: VecCoeficientes,VecRaices
real(8) tol

! COMIENZA EL PROGRAMA

CALL LeeCoeficientes(VecCoeficientes,Grado)
!----------ERROR--------!
write(*,*)"Ingrese la cota maxima de error" 
read(*,*)tol
!------Iteraciones-----!
write(*,*)"Ingrese maxima cantidad de iteraciones" 
read(*,*)maxiter


do while (CoeficienteCero(VecCoeficientes,Grado) .AND. Metodo<4) 

	Metodo=Metodo+1
	CALL PreProcesamiento(VecCoeficientes,Grado,Metodo,Diferencia) !hay 3 metodos para corregir 0, haciendo diferentes reempazos s=x+d s+cx s=1/x
  
enddo

if (Metodo>=4)

  write(*,*)"Imposible obtener raices del polinomio : No fue posible modificarlo para que sus coeficientes fueran distintos de 0")
  
 else
 
	VecRaices=0
	call QuotientDifferenceSub(VecCoeficientes,VecRaices,grado,tol,maxiter)
	call CorrigeSolucion(VecRaices,grado,Metodo,diferencia)
	
endif
!Write (*,*) 'Ingrese la tolerancia'               !!Aca esta en rojo porque mientras lo vas probando es más facil ya dejarlos puestos
                                                    ! están más abajo los valores.
   !Read(*,*) tol
   !error=2*tol
!Write (*,*) 'Ingrese la cantidad maxima de iteraciones'
   !Read(*,*) maxiter










CONTAINS

function CoeficienteCero(VecCoeficientes,Grado)
logical CoeficienteCero
real(8) VecCoeficientes(:)
integer Grado,i

do while ((VecCoeficientes(i)/=0) .AND. (i<Grado))
	i=i+1
enddo

CoeficienteCero=(Grado==i)



endfunction
SUBROUTINE LeeCoeficientes(VecCoeficientes,Grado)
integer i,Grado
real(8), allocatable :: VecCoeficientes(:)

open(2,FILE="CoeficientesPolinomio.txt", action="read")
read(2,*)grado

allocate(VecCoeficientes(grado))

do i=1,grado 
	read(2,*)VecCoeficientes(i)
enddo



END SUBROUTINE

subroutine QuotientDifferenceSub(VecCoeficientes,VecRaices,tol,maxiter)
!parametros 
real(8) intent(inout)::VecCoeficientes(:),VecRaices(:)
integer ,intent(in):: grado,maxiter
real(8) ,intent(in)::error
!var locales
integer i
real(8) error

   
!COMIENZA QD
q=0.          !que carajo significan las variables
e(n)=0.
e(0)=0.


!CALL Iteracion_cero(q,n,e,i) 


 q(1)=-(a(n-1))/(a(n))
    do i=1, n-1
      e(i)=a(n-i-1)/a(n-i)
    end do
Write(*,*)' '

maxiter=20
iter=0
tol=0.0005
error=2*tol



Do while (error>tol .and. maxiter>iter)

    Do i=1, n                  !Calculo la fila de entera en cada iteración
       q(i)=e(i)-e(i-1)+q(i)     
    End do
    
    Do i=1, n-1                !Como e(0)=e(n)=0 siempre, uso esos límites en el ciclo
       e(i)=e(i)*q(i+1)/q(i)   ! Hay un error de tipeo en el pdf que figura q(i-1) pero en la tabla utiliza q(i+1), sino no da
    end do
   

     iter=iter+1
     error=maxval(abs(e))
end do

end subroutine

End program
