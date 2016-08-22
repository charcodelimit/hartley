;;; -*- Mode: Lisp; Syntax: COMMON-LISP; Package: ; Base: 10 -*-
;;; -*- coding: utf-8 -*-
;;;****************************************************************************
;;; FILE:        fht.lisp
;;; LANGUAGE:    Common-Lisp
;;;
;;; DESCRIPTION
;;;  Fast Hartley Transformation implemented in Common Lisp
;;;
;;;
;;; Author: Christian Hofmann-Fuchs
;;;
;;; Created: Di Okt 12 11:29:21 2004 (+0200)
;;;
;;; Last-Updated: Mo Aug 22 22:14:53 2016 (+0200)
;;;           By: Christian Hofmann-Fuchs
;;;           Update #: 4
;;;
;;; Copyright (C) 2004,2016, Christian Hofmann-Fuchs. All rights reserved.
;;;
;;; Permission is hereby granted, free of charge, to any person
;;; obtaining a copy of this software and associated documentation
;;; files (the "Software"), to deal in the Software without
;;; restriction, including without limitation the rights to use,
;;; copy, modify, merge, publish, distribute, sublicense, and/or
;;; sell copies of the Software, and to permit persons to whom the
;;; Software is furnished to do so, subject to the following
;;; conditions:
;;;
;;; The above copyright notice and this permission notice shall be
;;; included in all copies or substantial portions of the Software.
;;;
;;; THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND,
;;; EXPRESS OR IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES
;;; OF MERCHANTABILITY, FITNESS FOR A PARTICULAR PURPOSE AND
;;; NONINFRINGEMENT. IN NO EVENT SHALL THE AUTHORS OR COPYRIGHT
;;; HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER LIABILITY,
;;; WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING
;;; FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR
;;; OTHER DEALINGS IN THE SOFTWARE.
;;;
;;;****************************************************************************

(in-package :cl-user)

(defclass transformable()

  ((data :type '(simple-array float *)
	 :accessor data
	 :initarg :data
	 :documentation "The storage containing the data to be transformed.")
   
   (accu :type '(simple-array float 2)
	 :accessor accu
	 :initarg :accu
	 :documentation "A storage for temporary transformation results.")

   (len :type fixnum
	:accessor len
	:initarg :len
	:documentation "The length of the data vector.")

   (sne  :type (simple-array float *)
	 :accessor sne
	 :initarg :sne
	 :documentation "A sine table for fast lookup.")

   (csn  :type (simple-array float *)
	 :accessor csn
	 :initarg :csn
	 :documentation "A cosine table for fast lookup."))

  (:default-initargs :data (make-array 0 :element-type 'float)
		     :len 0)
  
  (:documentation "A transformable contains data which can be transformed.")
  )

(defun new-transformable (simple-float-vector)
  (if (eq simple-float-vector ())
      (make-instance 'transformable)
    (make-instance 'transformable :data simple-float-vector)))


;;; some utilities
(defmacro swap (x y)
  `(let ((tmp-var ,x))
     (setf ,x ,y)
     (setf ,y tmp-var)))


;;; trigonometric function lookup table initialization

(defmacro sne-r (instance x)
  `(svref (sne ,instance) ,x))

(defmacro csn-r (instance x)
  `(svref (csn ,instance) ,x))

;; chr: array offset clean
(defmethod create-trig-table ((tr transformable) number-of-points)
  (let ((angle 0)
	(omega (/ (* 2 pi) number-of-points)))
    ;; init tables
    (setf (sne tr) (make-array number-of-points :element-type 'float :initial-element 0.0))
    (setf (csn tr) (make-array number-of-points :element-type 'float :initial-element 0.0))
    ;; fill tables
    (dotimes (cnt number-of-points)
      (setf (sne-r tr cnt) (sin angle))
      (setf (csn-r tr cnt) (cos angle))
      (incf angle omega)))
  (setf (len tr) number-of-points))


;; fast hartley transformation
(defmacro create-accu (instance length)
  `(setf (accu ,instance) (make-array (list 2 ,length) :element-type 'float :initial-element 0.0)))

(defmacro accu-r (instance x y)
  `(aref (accu ,instance) ,x ,y))

(defun permute (index log-length)
  (let ((idx index)
	(j 0)
	(s 0))
    (dotimes (cnt log-length j)
      (setf s (ash idx -1))
      (setf j (+ (ash j 1) (- idx (ash s 1))))
      (setf idx s))))

(defmacro add-sub-butterfly (instance index1 index2 from-buffer to-buffer)
  `(let ((tmp (accu-r ,instance ,from-buffer ,index1)))
     (declare (float tmp))
     (setf (accu-r ,instance ,to-buffer ,index1)
	   (+ tmp (accu-r ,instance ,from-buffer ,index2)))
     (setf (accu-r ,instance ,to-buffer ,index2)
	   (- tmp (accu-r ,instance ,from-buffer ,index2)))))

(defmacro butterfly (instance index1 index2 angle from-buffer)
  `(let ((res-list '()))
     (push (+ (* (accu-r ,instance ,from-buffer ,index1)
		 (csn-r ,instance ,angle))
	      (* (accu-r ,instance ,from-buffer ,index2)
		 (sne-r ,instance ,angle)))
	   res-list)
     (push (- (* (accu-r ,instance ,from-buffer ,index1)
		 (sne-r ,instance ,angle))
	      (* (accu-r ,instance ,from-buffer ,index2)
		 (csn-r ,instance ,angle)))
	   res-list)
     res-list))

(defmacro plus-minus (instance val index index-offset next-array last-array)
  `((lambda () (setf (accu-r ,instance ,next-array ,index)
		     (+ (accu-r ,instance ,last-array ,index) ,val))
      (setf (accu-r ,instance ,next-array (+ ,index ,index-offset))
	    (- (accu-r ,instance ,last-array ,index) ,val)))))

(defmethod initialize ((tr transformable) src-vec-length)
  (create-trig-table tr src-vec-length)
  (create-accu tr src-vec-length))

;;;   The fast Hartley transform.
;; chr: log-length is also called stages in some implementations
;; chr: array-displacement clean
;; chr: should throw an exception if (eq (data tr) ())
(defmethod hartley-transform ((tr transformable))
  (declare (optimize speed))
  (let ((last-array 1) 
	(next-array 0) 
	(src-vec-length (length (data tr)))
	(stages (round (/
			(log (length (data tr)))
			(log 2)))))
    (declare (type integer last-array next-array src-vec-length stages))

    ;;(format t "[Starting]~%")

    (if (not (= (len tr) src-vec-length))
	(initialize tr src-vec-length))

    ;;(format t "[Initializing]~%")

    ;; permute elements
    (dotimes (cnt src-vec-length) (setf
				   (aref (accu tr) next-array (permute cnt stages))
				   (svref (data tr) cnt)))
    ;;(format t "[Add Sub]~%")
    
    ;; addition and subtraction
    (do ((cnt 0 (+ cnt 2)))
	((>= cnt (- src-vec-length 1)))
      (declare (integer cnt))
      (add-sub-butterfly tr cnt (+ cnt 1) next-array next-array))

    ;;(format t "~A~%" (accu tr))
    
    (do ((cnt 0 (+ cnt 4)))
	((>= cnt src-vec-length))
      (declare (integer cnt))
      (add-sub-butterfly tr cnt (+ cnt 2) next-array last-array)
      (add-sub-butterfly tr (+ cnt 1) (+ cnt 3) next-array last-array)
      ;;(format t "~A~%" (accu tr))
    )

    ;;(format t "[Looping]~%")
    
    (let ((s 4)
	  (u (- stages 1)))
      (declare (integer s u))
      (do ((cnt 2 (+ cnt 1)))
	  ((>= cnt stages))
	(declare (integer cnt))
	(decf u)
	(let ((s2 (ash s 1))
	      (s0 (ash 1 (- u 1)))) ; raise one to 2^log2(N)-3 
	  (do ((q 0 (+ s2 q))) ; multiples of 8
	      ((>= q src-vec-length))
	    (declare (integer q))
	    (let ((index q)
		  (k (- (+ q s) 1)))
	      (declare (integer index k))

	      ;;(format t "[Add-Sub butterfly]~%")
	      
	      (add-sub-butterfly tr index (+ q s) last-array next-array) ; last -> next
	      (do ((j s0 (+ j s0)))
		  ((> j (ash src-vec-length -2)))
		(declare (integer j))

		(incf index)

		;;(format t "[Butterfly]~%")
		
		(let ((res (butterfly tr (+ index s) (+ k s) j last-array)))
		  (plus-minus tr (cadr res) index s next-array last-array)
		  (plus-minus tr (car res) k s next-array last-array))
		
		(decf k))))
	  (swap next-array last-array)
	  (setf s s2))))

    ;;(format t "[Finishing]~%")
    
    (let ((res (make-array src-vec-length :element-type 'float :initial-element 0.0)))
      (dotimes (cnt src-vec-length)
    	(setf (svref res cnt) (accu-r tr last-array cnt)))
      (new-transformable res))))    
  
;;; testing

(defun line (length)
  (let ((test-data (make-array length :element-type 'float :initial-element 0.0)))
    (dotimes (cnt length)
      (setf (svref test-data cnt) 1))
    (scale test-data (/ 1 length))))

(defun sine-wave (length)
  (let ((angle 0)
	(omega (* 2 (/ pi length)))
	(test-data (make-array length :element-type 'float :initial-element 0.0)))
    (dotimes (cnt length)
      (setf (svref test-data cnt) (sin angle))
      (incf angle omega))
    test-data))

(defun square-wave (length)
  (let ((test-data (make-array length :element-type 'float :initial-element 0.0)))
    (dotimes (cnt (/ length 2))
      (setf (svref test-data cnt) 1))
     (dotimes (cnt (/ length 2))
      (setf (svref test-data (+ cnt (/ length 2))) 0))
     test-data))

;;; some helpers

;; returns a scaled vector
(defun scale (vector factor)
  "scale vector factor

Multiplies the contents of vector with the value given by parameter factor
A newly created vector with the scaled values is returned."
  
  (let ((new-vector (make-array (length vector) :element-type 'float :initial-element 0.0)))
    (dotimes (cnt (length vector))
      (setf (svref new-vector cnt) (* factor (svref vector cnt))))
    new-vector))

(defun plot-spectrum (vector &key (title "Power Spectrum") (screen-width 64) (stream t))
  (let* ((maximum (reduce #'max vector))
         (minimum (reduce #'min vector))
         (digits (length (write-to-string (length vector))))
         (width (- screen-width (+ digits 2)))
         (value 0)
         (control-string (concatenate 'string "~%~" (write-to-string digits) ",'0d ")))
    ;; title
    (terpri)
    (dotimes (i screen-width) (format stream "="))
    (format stream "~% ~A~%" title)
    (dotimes (i screen-width) (format stream "="))
    ;; Spectrum
    (dotimes (i (length vector))
      (setq value
            (round (* (/ (aref vector i) (- maximum minimum)) width)))
      (format stream control-string i)
      (if (< value 0)
          (progn
            (dotimes (j (- value (round (* (/ minimum (- maximum minimum)) width))))
              (format t " "))
            (dotimes (j (abs value))
              (format t "*"))
            (format stream "|"))
          (progn
            (dotimes (j (round (* (/ (abs minimum) (- maximum minimum)) width)))
              (format stream " "))
            (format stream "|")
            (dotimes (j value)
              (format stream "*")))))
    (terpri)
    (dotimes (i screen-width) (format stream "_"))))

(defun test ()
  (let ((input (new-transformable (line 16))))
    (plot-spectrum (data (hartley-transform input))
                   :title "Power Spectrum - Line"))
  (let ((input (new-transformable (sine-wave 16))))
    (plot-spectrum (data (hartley-transform input))
                   :title "Power Spectrum - Sine Wave"))
  (let ((input (new-transformable (square-wave 16))))
    (plot-spectrum (data (hartley-transform input))
                   :title "Power Spectrum - Square Wave")))

;; CL-USER> (test)
;;
;; ================================================================
;;  Power Spectrum - Line
;; ================================================================
;; 00 |************************************************************
;; 01 |
;; 02 |
;; 03 |
;; 04 |
;; 05 |
;; 06 |
;; 07 |
;; 08 |
;; 09 |
;; 10 |
;; 11 |
;; 12 |
;; 13 |
;; 14 |
;; 15 |
;; ________________________________________________________________
;; ================================================================
;;  Power Spectrum - Sine Wave
;; ================================================================
;; 00                               |
;; 01                               |******************************
;; 02                               |
;; 03                               |
;; 04                               |
;; 05                               |
;; 06                               |
;; 07                               |
;; 08                               |
;; 09                               |
;; 10                               |
;; 11                               |
;; 12                               |
;; 13                               |
;; 14                               |
;; 15 ******************************|
;; ________________________________________________________________
;; ================================================================
;;  Power Spectrum - Square Wave
;; ================================================================
;; 00                     |****************************************
;; 01                     |******************************
;; 02                     |
;; 03                     |************
;; 04                     |
;; 05                     |********
;; 06                     |
;; 07                     |******
;; 08                     |
;; 09                     |****
;; 10                     |
;; 11                     |**
;; 12                     |
;; 13                   **|
;; 14                     |
;; 15 ********************|
;; ________________________________________________________________
