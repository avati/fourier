#!/usr/bin/guile -s
!#

(use-modules (ice-9 debugger))

(define pi 3.1415926535897932385)

(define (list->list-with-index vanilla-list)
  (let ((index -1))
    (map (lambda (x)
	   (set! index (+ index 1))
	   (list index x))
	 vanilla-list)))

(define (reduce func a-list seed)
  (cond ((null? a-list) seed)
	(else (reduce func 
		      (cdr a-list) 
		      (func (car a-list)
			    seed)))))

;; example - (DFT (list (cons x1 y1) (cons x2 y2)))
;; x1,y1 x2,y2 are real/imaginary pairs

(define (DFT list-in-time-domain)
  (let ((list-length (length list-in-time-domain))
	(indexed-list (list->list-with-index list-in-time-domain)))
    (map (lambda (i)
	   (let ((arg (* -1 2 pi (car i) (/ 1 list-length)))
		 (i-index (car i)))
	     (reduce +
		     (map (lambda (k)
			    (let* ((k-index (car k))
				   (cosarg 
				    (cos (* arg k-index)))
				   (sinarg 
				    (sin (* arg k-index)))
				   (x1 (real-part (cadr k)))
				   (y1 (imag-part (cadr k))))
				      (+
				       (/ (- (* x1 cosarg)
					     (* y1 sinarg))
					  list-length)
				       (* (sqrt -1)
					  (/ (+ (* x1 sinarg)
						(* y1 cosarg))
					     list-length)))))
			  indexed-list)
		     0)))
	 indexed-list)))

(define (DFT-inverse list-in-freq-domain)
  (let ((list-length (length list-in-freq-domain))
	(indexed-list (list->list-with-index list-in-freq-domain)))
    (map (lambda (i)
	   (let ((arg (*  2 pi (car i) (/ 1 list-length)))
		 (i-index (car i)))
	     (reduce +
		     (map (lambda (k)
			    (let* ((k-index (car k))
				   (cosarg 
				    (cos (* arg k-index)))
				   (sinarg 
				    (sin (* arg k-index)))
				   (x1 (real-part (cadr k)))
				   (y1 (imag-part (cadr k))))
			      (+
			       (- (* x1 cosarg)
				  (* y1 sinarg))
			       (* (sqrt -1)
				  (+ (* x1 sinarg)
				     (* y1 cosarg))))))
			  indexed-list)
		     0)))
	 indexed-list)))

(define range-cache-list '())
(define range-cache-vector '#())

(define (range-cache-list-build end num this-list)
  (cond ((= end num) #t)
	(else (begin
		(set! range-cache-list (append range-cache-list
					       (list this-list)))
		(range-cache-list-build end (+ num 1)
					(append this-list (list num)))))))

(range-cache-list-build 1025 0 '())
(set! range-cache-vector (list->vector range-cache-list))

(define (range-iter x ans)
    (cond ((= x 0) ans)
	  (else (range-iter (- x 1) (append (list (- x 1)) ans)))))

;(define (range x)
;   (range-iter x '()))

(define (range x)
  (vector-ref range-cache-vector x))

(define (bit-rev-iter num bitwidth ans)
  (cond ((= bitwidth 0)
	 ans)
	((= (modulo num 2) 0)
	 (bit-rev-iter (/ num 2)
		       (- bitwidth 1)
		       ans))
	(else
	 (bit-rev-iter (/ (- num 1) 2)
		       (- bitwidth 1)
		       (+ ans (expt 2 (- bitwidth 1)))))))

(define bit-rev-cache-list '())
(define bit-rev-cache-vect '#())

(define (bit-rev-cache-list-build end width)
  (cond ((= end 0) (set! bit-rev-cache-list
			 (append (list 0) bit-rev-cache-list)))
	(else
	 (set! bit-rev-cache-list (append (list (bit-rev-iter end width 0))
				   bit-rev-cache-list))
	 (bit-rev-cache-list-build (- end 1) width))))


(define (bit-rev num bitwidth)
  (bit-rev-iter num bitwidth 0))

(define (bit-rev num bitwidth)
  (vector-ref bit-rev-cache-vect num))

(define (log2-unsafe-iter num ans)
  (cond ((<= num 0) 1)
        ((= num 1) ans)
        (else (log2-unsafe-iter (/ num 2) (+ ans 1)))))

(define (log2-unsafe num)
  (log2-unsafe-iter num 0))

(define (log2 num)
  (/ (log num) (log 2)))

(define log2 log2-unsafe)

(define (FFT list-in-time-domain)
  (let* ((len (length list-in-time-domain))
	 (temp-vect (list->vector list-in-time-domain))
	 (vect (list->vector list-in-time-domain))
	 (base2 (log2 len))
	 (c1 -1)
	 (c2 0)
	 (l2 1))

    ;; bit reverse indexes
    (for-each (lambda (idx)
		(vector-set! vect (bit-rev idx base2)
			     (vector-ref temp-vect idx)))
	      (range len))

    (for-each 
     (lambda (l)
       (let* ((l1 l2)
	      (u1 1.0)
	      (u2 0.0))
	 (set! l2 (* 2 l2))
	 (for-each 
	  (lambda (j)
	    (for-each 
	     (lambda (i)
	       (let* ((i1 (+ i l1))
		      (xyi (vector-ref vect i))
		      (xyi1 (vector-ref vect i1))
		      (reali1 (real-part xyi1))
		      (imagi1 (imag-part xyi1))
		      (ti (+ (- (* u1 reali1)
				(* u2 imagi1))
			     (* (sqrt -1)
				(+ (* u1 imagi1)
				   (* u2 reali1))))))
		 (vector-set! vect i1 (- xyi ti))
		 (vector-set! vect i (+ xyi ti))))
	     (map (lambda (_num)
		    (+ (* _num l2) j))
		  (range (quotient (+ l2 -1 (- len j)) l2))))
	    (let ((z (- (* c1 u1) (* c2 u2))))
	      (set! u2 (+ (* u1 c2) (* u2 c1)))
	      (set! u1 z)))
	  (range l1))
	 (set! c2 (* -1 (sqrt (/ (- 1 c1) 2))))
	 (set! c1 (sqrt (/ (+ 1 c1) 2)))))
       (range base2))
  (map (lambda (n)
	 (/ n len))
       (vector->list vect))))

(define (FFT-inverse list-in-time-domain)
  (let* ((len (length list-in-time-domain))
	 (temp-vect (list->vector list-in-time-domain))
	 (vect (list->vector list-in-time-domain))
	 (base2 (log2 len))
	 (c1 -1)
	 (c2 0)
	 (l2 1))

    ;; bit reverse indexes
    (for-each (lambda (idx)
		(vector-set! vect (bit-rev idx base2)
			     (vector-ref temp-vect idx)))
	      (range len))

    (for-each 
     (lambda (l)
       (let* ((l1 l2)
	      (u1 1.0)
	      (u2 0.0))
	 (set! l2 (* 2 l2))
	 (for-each 
	  (lambda (j)
	    (for-each 
	     (lambda (i)
	       (let* ((i1 (+ i l1))
		      (xyi (vector-ref vect i))
		      (xyi1 (vector-ref vect i1))
		      (reali1 (real-part xyi1))
		      (imagi1 (imag-part xyi1))
		      (ti (+ (- (* u1 reali1)
				(* u2 imagi1))
			     (* (sqrt -1)
				(+ (* u1 imagi1)
				   (* u2 reali1))))))
		 (vector-set! vect i1 (- xyi ti))
		 (vector-set! vect i (+ xyi ti))))
	     (map (lambda (_num)
		    (+ (* _num l2) j))
		  (range (quotient (+ l2 -1 (- len j)) l2))))
	    (let ((z (- (* c1 u1) (* c2 u2))))
	      (set! u2 (+ (* u1 c2) (* u2 c1)))
	      (set! u1 z)))
	  (range l1))
	 (set! c2 (sqrt (/ (- 1 c1) 2)))
	 (set! c1 (sqrt (/ (+ 1 c1) 2)))))
       (range base2))
       (vector->list vect)))


(bit-rev-cache-list-build 1024 10)
(set! bit-rev-cache-vect (list->vector bit-rev-cache-list))


(for-each (lambda (idx)
	    (FFT-inverse (map (lambda (x) (* x 2)) (FFT (range 1024))))
	    (FFT-inverse (map (lambda (x) (* x 3)) (FFT (range 1024)))))
	  (range 100))
