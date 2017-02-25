(ns ssm4clj.core
  (:require [ssm4clj.misc :refer :all]
            [clojure.core.matrix :refer :all]
            [clojure.data.avl :as avl]))

(defn corr
 "Computes correlation matrix for length scale ls and delay d"
  [ls d]
  (let [l (/ (sqrt 5) ls)
        q (/ (* 16 (pow l 5)) 3)
        em2dl (exp (* l d -2))]
    (matrix [
             [(/ (* q (+ 3 (* em2dl (+ -3 (* -2 d l (+ 3 (* d l (+ 3 (* d l (+ 2 (* d l))))))))))) (* 16 (pow l 5)))
              (* 0.125 em2dl q (pow d 4))
              (* -1 (/ (* q (+ (* -1 em2dl) 1 (* 2 d l em2dl (+ -1 (* d l (+ -1 (* d l (+ -2 (* d l))))))))) (* 16 (pow l 3))))]
             [(* 0.125 em2dl q (pow d 4))
              (/ (* q (+ 1 (* em2dl (+ -1 (* -2 d l (+ 1 (* d l (pow (+ -1 (* d l)) 2)))))))) (* 16 (pow l 3)))
              (* 0.125 em2dl q (pow d 2) (pow (+ -2 (* d l)) 2))]
             [(* -1 (/ (* q (+ (* -1 em2dl) 1 (* 2 d l em2dl (+ -1 (* d l (+ -1 (* d l (+ -2 (* d l))))))))) (* 16 (pow l 3))))
              (* 0.125 em2dl q (pow d 2) (pow (+ -2 (* d l)) 2))
              (/ (* q (+ 3 (* em2dl (+ -3 (* -2 d l (+ -5 (* d l (+ 11 (* d l (+ -6 (* d l))))))))))) (* 16 l))]])))

(defn reverse-corr
 "Computes reverse-correlation matrix for length scale ls and delay d"
  [ls d]
  (let [l (/ (sqrt 5) ls)
        q (/ (* 16 (pow l 5)) 3)
        em2dl (exp (* l d -2))]
    (matrix [
             [(/ (* q (+ 3 (* em2dl (+ -3 (* -2 d l (+ 3 (* d l (+ 3 (* d l (+ 2 (* d l))))))))))) (* 16 (pow l 5)))
              (* -1 0.125 em2dl q (pow d 4))
              (* -1 (/ (* q (+ (* -1 em2dl) 1 (* 2 d l em2dl (+ -1 (* d l (+ -1 (* d l (+ -2 (* d l))))))))) (* 16 (pow l 3))))]
             [(* -1 0.125 em2dl q (pow d 4))
              (/ (* q (+ 1 (* em2dl (+ -1 (* -2 d l (+ 1 (* d l (pow (+ -1 (* d l)) 2)))))))) (* 16 (pow l 3)))
              (* -1 0.125 em2dl q (pow d 2) (pow (+ -2 (* d l)) 2))]
             [(* -1 (/ (* q (+ (* -1 em2dl) 1 (* 2 d l em2dl (+ -1 (* d l (+ -1 (* d l (+ -2 (* d l))))))))) (* 16 (pow l 3))))
              (* -1 0.125 em2dl q (pow d 2) (pow (+ -2 (* d l)) 2))
              (/ (* q (+ 3 (* em2dl (+ -3 (* -2 d l (+ -5 (* d l (+ 11 (* d l (+ -6 (* d l))))))))))) (* 16 l))]])))

(defn innovation
  "Computes covariance matrix for variance ρ² length-scale l and delay Δ"
  [ρ² l Δ] (mmul ρ² (corr l Δ)))

(defn approx-Qnvi
 "Returns unscaled Covariance matrix for length scale ls and delay d"
  [innovation ρ² ls d]
  (let [l (/ (sqrt 5) ls)
        q (/ (* 16 (pow l 5)) 3)]
    (if (< (* d l) 0.01 )
      (matrix [ [0 0 0] [0 0 0] [0 0 (/ 1 (* ρ² q d))]])
      (inverse (innovation d)))))

(defn reverse-innovation
  "Scales the reverse covariance matrix by variance parameter ρ²"
  [ρ² l Δ] (mmul ρ² (reverse-corr l Δ)))

(defn regression
  "Returns a regression matrix for length scale ls and delay d"
  [ls d]
  (let [l (/ (sqrt 5) ls)
        emdl (exp (* l d -1))]
    (matrix [
             [(* 0.5 emdl (+ 2 (* 2 d l) (* d d l l)))
              (* emdl d (+ 1 (* d l)))
              (* 0.5 emdl d d)]
             [(* -1 0.5 emdl d d l l l)
              (* -1 emdl (+ -1 (* -1 d l) (* d d l l)))
              (* -1 0.5 emdl d (+ -2 (* d l)))]
             [(* 0.5 emdl d l l l (+ -2 (* d l)))
              (* emdl d l l (+ -3 (* d l)))
              (* 0.5 emdl (+ 2 (* -4 d l) (* d d l l)))]])))

(defn reverse-regression
  "Returns a regression matrix for length scale ls and delay d"
  [ls d]
  (let [l (/ (sqrt 5) ls)
        emdl (exp (* l d -1))]
    (matrix [
             [(* 0.5 emdl (+ 2 (* 2 d l) (* d d l l)))
              (* -1 emdl d (+ 1 (* d l)))
              (* 0.5 emdl d d)]
             [(* 0.5 emdl d d l l l)
              (* -1 emdl (+ -1 (* -1 d l) (* d d l l)))
              (* 0.5 emdl d (+ -2 (* d l)))]
             [(* 0.5 emdl d l l l (+ -2 (* d l)))
              (* -1 emdl d l l (+ -3 (* d l)))
              (* 0.5 emdl (+ 2 (* -4 d l) (* d d l l)))]])))

(defn statcorr
  "Returns the stationary correlation matrix for length scale ls"
  [ls]
  (let [l (/ (sqrt 5) ls)
        q (/ (* 16 (pow l 5)) 3)]
      (matrix [
               [(/ (* 3 q) (* 16 (pow l 5)))
                0
                (* -1 (/ q (* 16 l l l)))]
               [0
                (/ q (* 16 l l l))
                0]
               [(* -1 (/ q (* 16 l l l)))
                0
                (/ (* 3 q) (* 16 l))]])))

(defn statcov
  "Computes the stationary covariance matrix for variance ρ²
  and length scale l"
  [ρ² l] (mmul ρ² (statcorr l)))

(defn sample-neighbour
 "Sample F[i+1] | F[i] or F[i-1] | F[i] depending on context"
 [F regression innovation delay]
 (let [X (regression delay)
       q (sample-safe-mvn (innovation delay))]
  (add (mmul X F) q)))

(defn sample-sandwich
 "Sample F[i] | F[i-1], F[i+1]"
 [Fv delay-v Fn delay-n regression innovation approx-Qnvi]
 (let [Xv (regression delay-v)
       Qv (innovation delay-v)
       Xn (regression delay-n)
       delay-nv (+ delay-n delay-v)
       Xnv (regression delay-nv)
       Qnvi (approx-Qnvi delay-nv)
       XvFv (mmul Xv Fv)
       XnvFv (mmul Xnv Fv)
       QvXn' (mmul Qv (transpose Xn))
       mean-vector (add XvFv (mmul QvXn' Qnvi (sub Fn XnvFv)))
       cov-matrix (sub Qv (mmul QvXn' Qnvi (transpose QvXn')))]
   (add mean-vector (sample-safe-mvn cov-matrix))))

(defn sample-gp
  "Sample a gaussian process function with specified variance
   and time-scale. Optionally provide known function values at
   a collection of times so that drawn function passes through
   these points."
  ([gp-var gp-time-scale]
   (sample-gp [] gp-var gp-time-scale))
  ([known-values gp-var gp-time-scale]
   (let [cond-set (atom (into (avl/sorted-map) known-values))
         reg (partial regression gp-time-scale)
         inn (partial innovation gp-var gp-time-scale)
         rev-reg (partial reverse-regression gp-time-scale)
         rev-inn (partial reverse-innovation gp-var gp-time-scale)
         app-Qnvi (partial approx-Qnvi inn gp-var gp-time-scale)]
     (with-meta
       (fn [t]
         (if (contains? @cond-set t) (@cond-set t)
             (let [vor (avl/nearest @cond-set < t)
                   nach  (avl/nearest @cond-set > t)
                   has-vor (not (nil? vor))
                   has-nach  (not (nil? nach))
                   has-both  (and has-vor has-nach)]
               (cond
                 has-both (let [[tv Fv] vor
                                [tn Fn] nach
                                delay-v (- t tv)
                                delay-n (- tn t)
                                Ft (sample-sandwich Fv delay-v Fn delay-n reg inn app-Qnvi)]
                            (do (swap! cond-set assoc t Ft) Ft))
                 has-vor (let [[tv Fv] vor
                               delay (- t tv)
                               Ft (sample-neighbour Fv reg inn delay)]
                           (do (swap! cond-set assoc t Ft) Ft))
                 has-nach (let [[tn Fn] nach
                                delay (- tn t)
                                Ft (sample-neighbour Fn rev-reg rev-inn delay)]
                            (do (swap! cond-set assoc t Ft) Ft))
                 :else (let [Ft (sample-safe-mvn (statcov gp-var gp-time-scale))]
                         (do (swap! cond-set assoc t Ft) Ft))))))
       {:var gp-var :time-scale gp-time-scale :data cond-set}))))

(defn AR-params
  "Yields time differences, covariance matrices, regression matrices
   and stationary covariance matricex determining the autoregressive
   structure on F"
  [times gp-var gp-time-scale]
  (let [delays (map - (rest times) times)
        Qs (map (partial innovation gp-var gp-time-scale) delays)
        Xs (map (partial regression gp-time-scale) delays)
        P (statcov gp-var gp-time-scale) ]
    [delays Qs Xs P]))

(defn filter-params
  "Augments the AR-paremeters with additional useful quantities"
  [times obs obs-var gp-var gp-time-scale]
  (let [[delays Qs Xs P] (AR-params times gp-var gp-time-scale)
        mar-obs-var-1 (+ obs-var (mget P 0 0))
        minit (div (mul (slice P 1 0) (first obs)) mar-obs-var-1)
        HP (slice P 0 0)
        Minit (sub P (div (outer-product HP HP) mar-obs-var-1))]
    [delays Qs Xs P HP mar-obs-var-1 minit Minit]))

(defn forward-filter
  "Returns new mean and covariance matrix along with other
   useful quantities"
 [inputs params]
 (let [[m M _ _ _] inputs
       [y obs-var X Q] params
       Xm (mmul X m)
       mar-obs-mean (mget Xm 0)
       XMXQ (add Q (mmul X M (transpose X)))
       XMXQH (slice XMXQ 0 0)
       HXMXQH (mget XMXQ 0 0)
       mar-obs-var (+ obs-var HXMXQH)
       residual (- y mar-obs-mean)
       m-out (add Xm (div (mmul XMXQH residual) mar-obs-var))
       M-out (sub XMXQ (div (outer-product XMXQH XMXQH) mar-obs-var))]
   [m-out M-out mar-obs-mean mar-obs-var XMXQH]))

(defn backward-compose
  "Returns the mean and covariance matrix defining the distribution
  of F_i conditional on F_i+1 and Y_1:i. Mean is a function of F_i+1."
 [m M X Q]
 (let [XMXQ (add Q (mmul X M (transpose X)))
       XMXQi (inverse XMXQ)
       MX' (mmul M (transpose X))
       XM (mmul X M)
       Xm (mmul X m)
       mean-function (fn [F] (add m (mmul MX' XMXQi (sub F Xm))))
       cov-matrix (sub M (mmul MX' XMXQi XM))]
   [mean-function cov-matrix]))

(defn backward-smooth
  "Returns the mode of F[i] | F[i+1] , Y[1], ...,Y[i]"
  [F distribution]
  (let [[mean-function cov-matrix] distribution]
    (mean-function F)))

(defn backward-sample
  "Samples a draw from F[i] | F[i+1], Y[1], ..., Y[i]"
 [F distribution]
 (let [[mean-function cov-matrix] distribution]
   (sample-safe-mvn (mean-function F) cov-matrix)))

(defn FFBC
 "Forward Filtering Backward Compose Algorithm"
 [times obs obs-var gp-var gp-time-scale]
 (let [[delays Qs Xs P HP mar-obs-var-1 minit Minit]
       (filter-params times obs obs-var gp-var gp-time-scale)
       FF (reductions forward-filter
                      [minit Minit 0 mar-obs-var-1 HP]
                      (map vector (rest obs) (repeat obs-var) Xs Qs))
       BCn [(first (last FF)) (second (last FF))]
       BC (map backward-compose (map first FF) (map second FF) Xs Qs)]
   [BC BCn]))

(defn FFBS
 "Forward Filtering Backward Sampling Algorithm"
 [times obs obs-var gp-var gp-time-scale]
 (let [[BC BCn] (FFBC times obs obs-var gp-var gp-time-scale)
       Fn (sample-safe-mvn (first BCn) (second BCn))
       Fs (reverse (reductions backward-sample Fn (reverse BC)))]
   (zipmap times Fs)))

(defn FFBM
  "Forward Filtering Backward Smoothing Algorithm"
  [times obs obs-var gp-var gp-time-scale]
  (let [[BC BCn] (FFBC times obs obs-var gp-var gp-time-scale)
        Fn (first BCn)
        Fs (reverse (reductions backward-smooth Fn (reverse BC)))]
    (zipmap times Fs)))

(defn logpdf-prior-gp
  "Evaluates the prior of the gp at F"
  ([F]
   (let [{gp-var :var
          gp-time-scale :time-scale} (meta F)]
     (logpdf-prior-gp F gp-var gp-time-scale)))
  ([F gp-var gp-time-scale]
   (let [data (deref (:data (meta F)))
         times (keys data)
         Fs (vals data)
         [delays Qs Xs P] (AR-params times gp-var gp-time-scale)
         logpdf-AR (fn [Fn Fv X Q] (logpdf-mvnormal Fn (mmul X Fv) Q))]
     (reduce +
             (logpdf-mvnormal (first Fs) (matrix [0 0 0]) P)
             (map logpdf-AR (rest Fs) Fs Xs Qs)))))

(defn logpdf-gen-F
 [F times obs obs-var gp-var gp-time-scale]
 (let [[BC BCn] (FFBC times obs obs-var gp-var gp-time-scale)
       Fs (map F times)
       Fn (last Fs)
       mean-functions (map first BC)
       means (map (fn [f x] (f x)) mean-functions (rest Fs))
       covars (map second BC)]
   (+ (reduce + (map logpdf-mvnormal (drop-last Fs) means covars))
      (logpdf-mvnormal Fn (first BCn) (second BCn)))))

(defn logpdf-fc-F
 [F times obs obs-var gp-var gp-time-scale]
 (let [[delays Qs Xs P HP mar-obs-var-1 minit Minit]
       (filter-params times obs obs-var gp-var gp-time-scale)
       Fs (map F times)
       means (cons [0 0 0] (map (fn [X F] (mmul X F)) Xs (drop-last Fs)))
       covars (cons P Qs)
       logdensity-f (reduce + (map logpdf-mvnormal Fs means covars))
       logdensity-y (reduce + (map logpdf-normal obs (map first Fs) (repeat obs-var)))]
   (+ logdensity-y logdensity-f)))

(defn log-likelihood
  "Evaluates the log likelihood having marginalized out the Gaussian Process
  with respect to its prior."
  [times obs obs-var gp-mean gp-var gp-time-scale]
  (if (or (neg? gp-var) (neg? gp-time-scale) (neg? obs-var))
    (Double/NEGATIVE_INFINITY)
    (let [[delays Qs Xs P HP mar-obs-var-1 minit Minit]
          (filter-params times obs obs-var gp-var gp-time-scale)
          FF (reductions forward-filter
                         [minit Minit 0 mar-obs-var-1 HP]
                         (map vector (rest obs) (repeat obs-var) Xs Qs))
          mar-means (map (fn [x] (nth x 2)) FF)
          mar-vars (map (fn [x] (nth x 3)) FF)]
      (reduce + (map (fn [x] (let [[y mu s2] x] (logpdf-normal y mu s2)))
                     (map vector obs mar-means mar-vars))))))

(defn mean-conditional
  "Computes the full condtional of the mean parameter with respect
  to a uniform prior having marginalized out the gaussian process
  function f with respect to its prior."
 [times obs obs-var gp-var gp-time-scale]
 (let
  [[delays Qs Xs P HP mar-obs-var-1 minit Minit]
   (filter-params times obs obs-var gp-var gp-time-scale)
   Cinit (mul (slice P 1 0) (/ (first obs) mar-obs-var-1))
   Dinit (negate (div (slice P 1 0) mar-obs-var-1))
   forward-filter
    (fn [acc params]
     (let
      [[pm prec C D m M] acc
       [y X Q] params
       [m+1 M+1 _ mar-obs-var XMXQH] (forward-filter [m M 0 0 0] [y obs-var X Q])
       XD (mmul X D)
       HXD (mget XD 0)
       HXD+1 (inc HXD)
       XC (mmul X C)
       HXC (mget XC 0)
       prec+1 (+ prec (/ (square HXD+1) mar-obs-var))
       pm+1 (+ pm (/ (* (sub y HXC) HXD+1) mar-obs-var))
       F (div XMXQH mar-obs-var)
       C+1 (add XC (negate (mmul F HXC)) (mmul F y))
       D+1 (sub XD F (mmul F HXD))]
      [pm+1 prec+1 C+1 D+1 m+1 M+1]))

   ffs (reduce forward-filter
               [(/ (first obs) mar-obs-var-1) (/ 1 mar-obs-var-1) Cinit Dinit minit Minit]
               (map vector (rest obs) Xs Qs))
   precision (second ffs)
   precxmean (first ffs)]
  [(/ precxmean precision) (/ 1.0 precision)]))
