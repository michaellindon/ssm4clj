(ns ssm4clj.stability
  (:require [ssm4clj.core :refer :all]
            [incanter.charts :as ip]
            [Incanter.core :as ic]
            [clojure.core.matrix.random :refer :all]
            [clojure.core.matrix :refer :all]))

(set-current-implementation :vectorz)

;;Generate a random stream of uniforms to work with
(def seed 1)
(def random-stream (randoms seed))
(def F (sample-gp [] 1 1))
(ic/view (ip/scatter-plot (take 1000 random-stream) (map (comp first F) (take 1000 random-stream))))

(ic/view (ip/scatter-plot (take 2000 random-stream) (map (comp first F) (take 2000 random-stream))))

(ic/view (ip/scatter-plot (take 5000 random-stream) (map (comp first F) (take 5000 random-stream))))

(ic/view (ip/scatter-plot (take 10000 random-stream) (map (comp first F) (take 10000 random-stream))))

(ic/view (ip/scatter-plot (take 100000 random-stream) (map (comp first F) (take 100000 random-stream))))

(def times (sort (take 100000 random-stream)))
(def obs (map (fn [x] (+ x (first (sample-normal 1)))) (map (comp first F) times)))
(ic/view (ip/scatter-plot times obs))
(def G (sample-gp (FFBS times obs 1 1 1) 1 1))

(ic/view (ip/scatter-plot times (map (comp first G) times)))
