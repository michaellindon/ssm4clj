(ns ssm4clj.misc
  (:require [clojure.core.matrix :refer :all]
            [clojure.core.matrix.linear :as la]
            [incanter.stats :refer :all]))

(defn logpdf-mvnormal
  "log density of multivariate normal"
  ([observations cov-matrix]
   (let [n (ecount observations)]
     (logpdf-mvnormal (zero-array [n]) cov-matrix)))
  ([observations mean-vector cov-matrix]
   (let [n (ecount observations)
         residual (sub observations mean-vector)
         exponent (* 0.5 (dot residual (mmul (inverse cov-matrix) residual)))
         normalizer (+ (* 0.5 (log (det cov-matrix)))
                       (* n 0.5 (log (* 2 Math/PI))))]
     (negate (+ normalizer exponent)))))

(defn logpdf-normal
  "log density of univariate normal"
 [y mu s2]
 (let [exponent (* 0.5 (/ (square (- y mu)) s2))
       normalizer (* 0.5 (log (* 2 Math/PI s2)))]
   (negate (+ exponent normalizer))))

(defn sample-safe-mvn
 "Numerically stable function to sample a multivariate normal"
  ([cov-matrix]
   (let [ncols (column-count cov-matrix)
         z (sample-normal ncols)
         d (la/svd cov-matrix)
         {U :U S :S V* :V*} (la/svd cov-matrix)
         D (diagonal-matrix (map (fn [x] (sqrt (max x 0))) S))]
     (mmul U D z)))
  ([mean-vector cov-matrix]
   (add mean-vector (sample-safe-mvn cov-matrix))))

(defn sample-truncated-normal [mu s2 a b]
  "Generate a univariate random normal truncated to interval [a,b]"
 (let [sd (sqrt s2)
       normalizer (- (cdf-normal b :mean mu :sd sd) (cdf-normal a :mean mu :sd sd))
       uniform (sample-beta 1)]
  (quantile-normal (+ (* uniform normalizer) (cdf-normal a :mean mu :sd sd))
                   :mean mu :sd sd)))

(defn sample-left-normal [mu s2]
  "Generate a univariate random normal truncated to (-inf,0)"
 (let [sd (sqrt s2)
       normalizer (cdf-normal 0 :mean mu :sd sd)
       uniform (sample-beta 1)]
  (quantile-normal (* uniform normalizer)
                   :mean mu :sd sd)))

(defn sample-right-normal [mu s2]
  "Generate a univariate random normal truncated to (0,Inf)"
 (let [sd (sqrt s2)
       normalizer (- 1 (cdf-normal 0 :mean mu :sd sd))
       uniform (sample-beta 1)]
  (quantile-normal (+ (* uniform normalizer) (cdf-normal 0 :mean mu :sd sd))
                   :mean mu :sd sd)))

(defn sample-slice
 "Sample a univariate unnormalized log-density g from position x with length L"
 [g w x]
 (let [y (+ (g x) (negate (sample-exp 1)))
       u (sample-beta 1)
       lower-bound (first (filter (fn [x] (< (g x) y))
                                  (iterate  (fn [x] (- x w))
                                            (- x (* u w)))))
       upper-bound (first (filter (fn [x] (< (g x) y))
                                  (iterate (fn [x] (+ x w))
                                           (+ x (* (- 1 u) w)))))]
   (loop [l lower-bound
          u upper-bound]
     (let [z (first (sample-uniform 1 :min l :max u))]
       (cond
         (> (g z) y) z
         (> z x) (recur l z)
         :else (recur z u))))))

