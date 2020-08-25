breed [atoms atom]
breed [fl-ends fl-end]

atoms-own [
  fx     ; x-component of force vector
  fy     ; y-component of force vector
  vx     ; x-component of velocity vector
  vy     ; y-component of velocity vector
  mass
  pinned?
]

globals [
  eps ; used in LJ force
  sigma ; used in LJ force
  cutoff-dist ; each atom is influenced by its neighbors within this distance
  dt ; time step for the velocity verlet algorithm
  sqrt-2-kb-over-m  ; constant
  cone-check-dist ; each atom links with neighbors within this distance
  prev-lattice-view ; the lattice view in the previous time step
  upper-left-fl ; upper left force line -- shear
  left-fl ; tension
  right-fl ; tension
  orig-length
  prev-length
  f-app-auto
  f-app-prev
  total-external-force
  reported-ex-force
  median-ycor
  top-neck-atoms
  bottom-neck-atoms
]

;;;;;;;;;;;;;;;;;;;;;;
;; Setup Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all
  set eps .07
  set sigma .907
  set cutoff-dist 5
  set dt .1 ;.05
  set sqrt-2-kb-over-m (1 / 20)
  set cone-check-dist 1.5
  setup-atoms-and-links
  init-velocity
  if lattice-view = "hide-atoms" [
    ask atoms [ hide-turtle ]
  ]
  reset-timer
  reset-ticks
end

to setup-atoms-and-links
  if force-mode = "Tension" and atoms-per-column mod 2 = 0 [ set atoms-per-column atoms-per-column + 1]
  create-atoms atoms-per-row * atoms-per-column [
    set shape "circle"
    set color blue
    set mass 1
    ifelse lattice-view = "small-atoms" [
      set size .6 ]
      [set size .9]
  ]
  let x-dist 1 ; the distance between atoms in the x direction
  let y-dist sqrt (x-dist ^ 2 - (x-dist / 2) ^ 2) ; the distance between rows
  let ypos (- atoms-per-column * y-dist / 2) ;the y position of the first atom
  let xpos (- atoms-per-row * x-dist / 2) ;the x position of the first atom
  let rnum 0 ; row number, starts at 0 for easy modulo division
  ask atoms [ ; setting up the HCP structure
    if xpos >= (atoms-per-row * x-dist / 2)  [ ; condition for starting new row
      set rnum rnum + 1
      set xpos (- atoms-per-row * x-dist / 2) + (rnum mod 2) * x-dist / 2
      set ypos ypos + y-dist
    ]
    setxy xpos ypos
    set xpos xpos + x-dist
  ]

  ; values used in assigning atom positions
  let ymax max [ycor] of atoms
  let xmax max [xcor] of atoms
  let ymin min [ycor] of atoms
  let xmin min [xcor] of atoms
  let median-xcor (median [xcor] of atoms)
  set median-ycor (median [ycor] of atoms)

  (ifelse force-mode = "Shear"[
    ask atoms with [ ((xcor = xmin or xcor = xmin + (1 / 2) or xcor = xmax or xcor = xmax - (1 / 2)) and ( ycor < median-ycor))] [
      set pinned? 1
    ]
    ]
    force-mode = "Tension"[ ; refine shape once you decide what ratio of neck to shoulder you want
      ask atoms with [xcor = min [xcor] of atoms] [die]
      set xmin min [xcor] of atoms
      ;ask atoms with [(ycor >= max [ycor] of atoms - 2 or ycor <= min [ycor] of atoms + 2) and xcor <= max [xcor] of atoms - 4 and xcor >= min [xcor] of atoms + 4] [die]
      ask atoms with [(ycor >= max [ycor] of atoms - 1 or ycor <= min [ycor] of atoms + 1) and xcor <= max [xcor] of atoms - 3.5 and xcor >= min [xcor] of atoms + 3.5] [die]
      ask atoms with [xcor = max [xcor] of atoms or xcor = max [xcor] of atoms - .5 ] [set pinned? 1]
      let max-ycor-neck max [ycor] of (atoms with [ xcor <= max [xcor] of atoms - 3.5 and xcor >= min [xcor] of atoms + 3.5])
      let min-ycor-neck min [ycor] of (atoms with [ xcor <= max [xcor] of atoms - 3.5 and xcor >= min [xcor] of atoms + 3.5])
      set top-neck-atoms atoms with [xcor <= max [xcor] of atoms - 3.5 and xcor >= min [xcor] of atoms + 3.5 and ycor = max-ycor-neck]
      set bottom-neck-atoms atoms with [xcor <= max [xcor] of atoms - 3.5 and xcor >= min [xcor] of atoms + 3.5 and ycor = min-ycor-neck]
    ]
  force-mode = "Compression" [
      ask atoms with [xcor = max [xcor] of atoms or xcor = max [xcor] of atoms - .5 ] [set pinned? 1]
    ]
  )

  if create-dislocation? [ ; creating the dislocation
    let curr-y-cor median [ycor] of atoms
    let curr-x-cor median [xcor] of atoms
    let ii 0
    while [ curr-y-cor <= ceiling (ymax) ] [
      ask atoms with [
        ycor <= curr-y-cor + y-dist * .75
        and ycor >= curr-y-cor ] [
        (ifelse xcor <= curr-x-cor + x-dist * .75
          and xcor >= curr-x-cor [ die ]
          xcor >= curr-x-cor [ set xcor xcor - .05 * ii ]
          xcor < curr-x-cor [ set xcor xcor + .05 * ii ])
        ]
      set curr-y-cor curr-y-cor + y-dist
      set curr-x-cor curr-x-cor - x-dist / 2
      set ii ii + 1
    ]
   ]

  ask atoms [
    let in-radius-atoms (other atoms in-radius cutoff-dist)
    update-links in-radius-atoms
  ]
  ask links [ ; stylizing/coloring links
    set thickness .25
    color-links
  ]

  ifelse force-mode != "Shear" [ ; is there a better way to do this fl assignment?
    create-fl-ends 2
    set left-fl xmin
    set right-fl xmax
    set orig-length right-fl - left-fl
    ask one-of fl-ends with [xcor = 0 and ycor = 0] [
      set xcor left-fl
      set ycor ymax + 2 ]
    ask one-of fl-ends with [xcor = 0 and ycor = 0] [
      set xcor left-fl
      set ycor ymin - 2
      create-link-with one-of other fl-ends with [xcor = left-fl]]
    ifelse force-mode = "Tension" [
      ask fl-ends [
        set color white
        set heading 270
      ]
    ]
    [
      ask fl-ends [
        set color white
        set heading 90
      ]
    ]
    set prev-length orig-length
    set f-app-auto f-app
  ]
  [
    create-fl-ends 2
    set upper-left-fl min [xcor] of atoms with [ ycor >= median [ycor] of atoms ]
    ask one-of fl-ends with [xcor = 0 and ycor = 0] [
      set xcor upper-left-fl
      set ycor ymax + 2 ]
    ask one-of fl-ends with [xcor = 0 and ycor = 0] [
      set xcor upper-left-fl
      set ycor median [ycor] of atoms
      hide-turtle
      create-link-with one-of other fl-ends]
    ask fl-ends [
      set color white
      set heading 90 ]
  ]
  ask links with [ is-fl-end? one-of both-ends ] [
      set color white
      ;set thickness .25
    ]
  ask atoms with [pinned? = 1] [ set shape "circle 2"]
end

to init-velocity ; initializes velocity for each atom based on the initial system-temp. Creates a random aspect in the
                 ; velocity split between the x velocity and the y velocity
  let speed-avg sqrt-2-kb-over-m * sqrt system-temp
  ask atoms [
    let x-portion random-float 1
    set vx speed-avg * x-portion * positive-or-negative
    set vy speed-avg * (1 - x-portion) * positive-or-negative]
end

to-report positive-or-negative
  report ifelse-value random 2 = 0 [-1] [1]
end


;;;;;;;;;;;;;;;;;;;;;;;;
;; Runtime Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;;;

to go
  if lattice-view != prev-lattice-view [ update-lattice-view ]
  set total-external-force 0
  control-temp
  ask links with [ is-atom? one-of both-ends ] [die]
  ask atoms [ ; moving happens before velocity and force update in accordance with velocity verlet
    move
  ]
  calculate-fl-positions
  if force-mode = "Tension" [ check-eq-adj-force ]
  ask atoms [
    update-force-and-velocity-and-links
  ]
  ask links with [is-atom? one-of both-ends ] [ ; stylizing/coloring links
    set thickness .25
    color-links
  ]
  set reported-ex-force total-external-force
  set prev-lattice-view lattice-view
  tick-advance dt
  update-plots
end

to update-lattice-view
  (ifelse lattice-view = "large-atoms" [
    ask atoms [
      show-turtle
      set size .9
    ]
  ]
  lattice-view = "small-atoms" [
    ask atoms [
       show-turtle
       set size .6
    ]
  ]
  [; lattice-view = hide-atoms
      ask atoms [ hide-turtle ]
  ])
end

to control-temp ; this heats or cools the system based on the average temperature of the system compared to the set system-temp
  let current-speed-avg mean [ sqrt (vx ^ 2 + vy ^ 2)] of atoms
  let target-speed-avg sqrt-2-kb-over-m * sqrt system-temp
  let scaling-factor target-speed-avg / current-speed-avg
  if current-speed-avg != 0 [
    ask atoms [
      set vx vx * scaling-factor
      set vy vy * scaling-factor
    ]
  ]
end

to move  ; atom procedure, uses velocity-verlet algorithm
  ifelse force-mode = "Shear" [
    if pinned? = 0 [
      set xcor velocity-verlet-pos xcor vx (fx / mass)
      set ycor velocity-verlet-pos ycor vy (fy / mass)]
    if xcor > max-pxcor [ ; kills force-arrows when their associated atoms move off the world
      die ; kills atoms when they move off the world
    ]
  ]
  [ ;; force-mode = "Tension" or "Compression"
    if pinned? = 0 [
      set xcor velocity-verlet-pos xcor vx (fx / mass)
      set ycor velocity-verlet-pos ycor vy (fy / mass)]
      if xcor > max-pxcor or xcor < min-pxcor [
        die ; kills atoms when they move off the world
      ]
   ]
end

to calculate-fl-positions
  ifelse force-mode = "Shear" [
    set upper-left-fl min [xcor] of atoms with [ ycor >= median-ycor ] ;+ .25 ]
    ask fl-ends [ set xcor upper-left-fl]
  ]
  [ ; force-mode = tension or compression
    set left-fl min [xcor] of atoms
    ask fl-ends with [xcor < 0] [ set xcor left-fl]
    ]
end

to check-eq-adj-force
  if f-app != f-app-prev [ set f-app-auto f-app ]
  if precision prev-length 4 = precision (right-fl - left-fl) 4 [ set f-app-auto f-app-auto + .01 ]
  set prev-length (right-fl - left-fl)
  set f-app-prev f-app
end

to update-force-and-velocity-and-links
  let new-fx 0
  let new-fy 0
  let total-force 0
  let in-radius-atoms (other atoms in-radius cutoff-dist)
  ask in-radius-atoms [ ; each atom calculates the force it feels from its neighboring atoms and sums these forces
    let r distance myself
    let force LJ-force r
    set total-force total-force + abs(force) ; keep track of this for visualization (coloring the atoms based on their potential energy -- since we do abs(force), that corresponds to PE)
    face myself
    set new-fx new-fx + (force * dx)
    set new-fy new-fy + (force * dy)
    ]

  ; adjusting the forces to account for any external applied forces
  let ex-force report-new-force
  set new-fx ex-force + new-fx
  set total-external-force total-external-force + ex-force

  ; updating velocity and force
  set vx velocity-verlet-velocity vx (fx / mass) (new-fx / mass)
  set vy velocity-verlet-velocity vy (fy / mass) (new-fy / mass)
  set fx new-fx
  set fy new-fy

  ifelse ex-force = 0 [update-atom-color total-force]
  [set color white]
  update-links in-radius-atoms
end

to update-atom-color [total-force] ; updating atom color
  (ifelse update-color? [
    set-color total-force
  ]
   [ set color blue ])
end

to update-links [in-radius-atoms] ; updating links
  if diagonal-right-links [
    set heading 330
    link-with-atoms-in-cone in-radius-atoms
  ]
  if diagonal-left-links [
    set heading 30
    link-with-atoms-in-cone in-radius-atoms
  ]
  if horizontal-links [
    set heading 90
    link-with-atoms-in-cone in-radius-atoms
  ]
end

to link-with-atoms-in-cone [atom-set]
  let in-cone-atoms (atom-set in-cone cone-check-dist 60)
    if any? in-cone-atoms [
      create-link-with min-one-of in-cone-atoms [distance myself]
    ]
end

to-report report-new-force ; change to external force only
  (ifelse force-mode = "Shear" [
    ifelse ycor >= median-ycor and xcor >= upper-left-fl and xcor <= upper-left-fl + .85 [
      ;report f-app * 1 / distancexy upper-left-fl ycor
       report f-app
    ]
    [ report 0]
    ]
  force-mode = "Tension" [
      let dist-l distancexy left-fl ycor
      ifelse dist-l <= 3.5 [ ; round vs no round?
       ; report -1 * f-app-auto * 1 / (dist-l + .5) + f-app-auto * 1 / (force-cutoff + .5)]
        ;report -1 * f-app-auto * 1 / round (dist-l + .51) + f-app-auto * 1 / round (force-cutoff + .51)
        report -1 * f-app-auto / 10]
        [report 0 ]
    ]
  force-mode = "Compression" [ ; this vs just force?
    let dist-l distancexy left-fl ycor
     ifelse dist-l <= 1 [ ; ceiling vs no ceiling?
        report (f-app * 1 / (dist-l + .5) - f-app * 1 / (1 + .5))]
        ;report f-app-auto * 1 / round (dist-l + .51) - f-app-auto * 1 / round (force-cutoff + .51)]
        ;report f-app-auto]
        [report 0 ]
    ]
  )
end

to-report LJ-force [ r ] ; optimize? + = attract, - = repulse (this derivative would usually have a negative sign in front, so it's as if we multiplied it by a negative)
  report (48 * eps / r )* ((sigma / r) ^ 12 - (1 / 2) * (sigma / r) ^ 6) + .0001
end

to-report velocity-verlet-pos [pos v a]  ; position, velocity and acceleration
  report pos + v * dt + (1 / 2) * a * (dt ^ 2)
end

to-report velocity-verlet-velocity [v a new-a]  ; velocity, acceleration, new acceleration
  report v + (1 / 2) * (new-a + a) * dt
end

to set-color [v]
  set color scale-color blue sqrt(v) -.3 1.7
end

to-report strain ; dc? true strain?
  report ((right-fl - left-fl) - orig-length) / orig-length
end

to-report stress ;!= 0.....ok?
;  let xcor-vals (range left-fl right-fl)
;  let min-A 10000000
;  let min-A-xcor 0
;  foreach xcor-vals [ x ->
;    let tmp-min-y [ycor] of min-one-of atoms with [xcor >= x - .5 and xcor <= x + .5] [ycor]
;    let tmp-max-y [ycor] of max-one-of atoms with [xcor >= x - .5 and xcor <= x + .5] [ycor]
;    if tmp-max-y - tmp-min-y < min-A and tmp-max-y - tmp-min-y != 0 [ set min-A-xcor x
;      set min-A tmp-max-y - tmp-min-y ]
;  ]
;  show min-A-xcor
;  ask atoms [ ifelse xcor >= min-A-xcor - .5 and xcor <= min-A-xcor + .5 [ set color red] [ set color blue]]
;  report (abs(total-external-force) / min-A) * 100

  let avg-max mean [ycor] of top-neck-atoms
  let avg-min mean [ycor] of bottom-neck-atoms
  let min-A avg-max - avg-min
  report (abs(total-external-force) / min-A) * 100
end

to color-links ; difficult to see......
  let min-eq-bond-len .995
  let max-eq-bond-len 1.018073
  (ifelse
    link-length < min-eq-bond-len [
      let tmp extract-rgb scale-color red sqrt(min-eq-bond-len - link-length) 1 -.2
      set color insert-item 3 tmp (125 + (1 + sqrt(min-eq-bond-len - link-length)) * 30) ]
    link-length > max-eq-bond-len [
      let tmp extract-rgb scale-color yellow sqrt (link-length - max-eq-bond-len) 1 -.2
      set color insert-item 3 tmp (125 + (1 + sqrt(link-length - max-eq-bond-len)) * 30)]
    [ let tmp extract-rgb white
      set color insert-item 3 tmp 125 ])
end

;to color-links ; pending update
;  let min-eq-bond-len .995
;  let max-eq-bond-len 1.018073
;  (ifelse
;    link-length < min-eq-bond-len [ set color scale-color red sqrt (min-eq-bond-len - link-length) -.05 .35 ]
;    link-length > max-eq-bond-len [ set color scale-color yellow sqrt (link-length - max-eq-bond-len) -.05 .35 ]
;    [ set color gray ])
;end

;to color-links ; difficult to see......
;  let min-eq-bond-len .995
;  let max-eq-bond-len 1.018073
;  (ifelse
;    link-length < min-eq-bond-len [
;      let tmp extract-rgb scale-color red sqrt(min-eq-bond-len - link-length) .9 0
;      set color insert-item 3 tmp (125 + (min-eq-bond-len - link-length) * 250) ]
;    link-length > max-eq-bond-len [
;      let tmp extract-rgb scale-color yellow sqrt (link-length - max-eq-bond-len) .9 0
;      set color insert-item 3 tmp (125 + (link-length - max-eq-bond-len) * 250)]
;    [ let tmp extract-rgb white
;      set color insert-item 3 tmp 125 ])
;end
@#$#@#$#@
GRAPHICS-WINDOW
192
10
820
639
-1
-1
20.0
1
10
1
1
1
0
1
0
1
-15
15
-15
15
1
1
1
ticks
30.0

BUTTON
11
194
97
227
NIL
setup
NIL
1
T
OBSERVER
NIL
NIL
NIL
NIL
1

BUTTON
105
194
190
227
NIL
go
T
1
T
OBSERVER
NIL
NIL
NIL
NIL
0

CHOOSER
25
13
163
58
force-mode
force-mode
"Shear" "Tension" "Compression"
0

SLIDER
12
288
184
321
system-temp
system-temp
0
.75
0.23
.01
1
NIL
HORIZONTAL

SLIDER
12
331
184
364
f-app
f-app
0
2
0.64
.01
1
N
HORIZONTAL

SWITCH
19
65
178
98
create-dislocation?
create-dislocation?
0
1
-1000

SWITCH
832
20
964
53
update-color?
update-color?
1
1
-1000

CHOOSER
14
236
181
281
lattice-view
lattice-view
"large-atoms" "small-atoms" "hide-atoms"
1

SWITCH
835
62
994
95
diagonal-right-links
diagonal-right-links
0
1
-1000

SWITCH
835
102
987
135
diagonal-left-links
diagonal-left-links
0
1
-1000

SWITCH
834
140
972
173
horizontal-links
horizontal-links
0
1
-1000

SLIDER
13
106
185
139
atoms-per-row
atoms-per-row
5
20
12.0
1
1
NIL
HORIZONTAL

SLIDER
13
148
185
181
atoms-per-column
atoms-per-column
5
20
10.0
1
1
NIL
HORIZONTAL

MONITOR
834
191
986
236
total external force
reported-ex-force
3
1
11

PLOT
834
254
1034
404
Stress-Strain Curve - Tension
strain
stress
0.0
0.1
0.0
20.0
true
false
"" ""
PENS
"default" 1.0 0 -16777216 true "" "if force-mode = \"Tension\" [ plotxy strain stress ]"

MONITOR
842
433
915
478
NIL
f-app-auto
5
1
11

MONITOR
958
435
1015
480
length
right-fl - left-fl
5
1
11

@#$#@#$#@
## WHAT IS IT?

This model displays the common natural phenomenon expressed by the Coulomb's inverse-square law.  It shows what happens when the strength of the force between two charges varies inversely with the square of the distance between them.

## HOW IT WORKS

In this model the formula used to guide each charge's behavior is the standard formula for Coulomb's law:

F = (q1 * q2 * Permittivity) / (r^2)

In this formula:
- "F" is the force between the two charges q1 and q2.
- "Permittivity", the constant of proportionality here, is a property of the medium in which the two charges q1 and q2 are situated.
- "r" is the distance between the centers of the two charges.

This is a single force two body model, where we have a charge q1 (the particle that is created when you press SETUP) and a proton (q2) (the blue particle that appears when you press the mouse button in the view).  If a particle is positively charged, it is colored blue. If it's negatively charged, it will be orange. The force is entirely one-way: only q1 is attracted towards (or repelled from) the proton (q2), while the proton (q2) remains unaffected.  Note that this is purely for purposes of simulation.  In the real world, Coulomb's force acts on all bodies around it.

Gravity is another example of an inverse square force.  Roughly speaking, our solar system resembles a nucleus (sun) with electrons (planets) orbiting around it.

For certain values of q1 (which you can control by using the CHARGE slider), you can watch q1 form elliptic orbits around the mouse pointer (q2), or watch q1 slingshot around q2, similar to how a comet streaks past our sun. The charges q1 and q2 are always set equal in magnitude to each other, although they can differ in their sign.

## HOW TO USE IT

When you press the SETUP button, the charge q1 is created in a medium determined by the permittivity value from the PERMITTIVITY slider. When you click and hold down the mouse anywhere within the view, the model creates a unit of positive charge (q2) at the position of the mouse.

The CHARGE slider sets the value of the charge on q1.  First, select the value of CHARGE on q1. You will see that the color of q1 reflects its charge. (For simulation ease, value of the charge on q2 is set to be the absolute value of this charge. Thus, it also determines at what distances the particles can safely orbit before they get sucked in by an overwhelming force.)

The FADE-RATE slider controls how fast the paths marked by the particles fade.  At 100% there won't be any paths as they fade immediately, and at 0% the paths won't fade at all.

The PERMITTIVITY slider allows you to change values of the constant of proportionality in Coulomb's law. What does this variable manipulate? The charges or the medium in which the charges are immersed?

When the sliders have been set to desirable levels, press the GO button to begin the simulation.  Move the mouse to where you wish q2 to begin, and click and hold the mouse button. This will start the particles moving. If you wish to stop the simulation (say, to change the value of CHARGE), release the mouse button and the particles will stop moving. You may then change any settings you wish. Then, to continue the simulation, simply put your mouse in the window again and click and hold. Objects in the window will only move while the mouse button is pressed down within the window.

## THINGS TO NOTICE

The most important thing to observe is the behavior of q1, the particle first placed in the world at SETUP.

What is the initial velocity for q1?

What happens as you change the value of q1 from negative to positive?

As you run the model, watch the graphs on the right hand side of the world. What can you infer from the graphs about the relationship between potential energy and distance between charges? What can you say about the relationship between Coulomb's force and distance between the charges from the graphs?

Move the mouse around and watch what happens if you move it quickly or slowly. Jiggle it around in a single place, or let it sit still. Observe what patterns the particles fall into. (You may keep FADE-RATE low to watch this explicitly.)

## THINGS TO TRY

Run the simulation playing with different values of:
a) charge - make sure to watch how different values of the CHARGE slider impact the model for any fixed value of permittivity.
b) permittivity - make sure to watch how different values of the PERMITTIVITY slider impact the model for any fixed value of charge.

Can you make q1 revolve around q2?  Imagine, if q1 would be an electron and q2 a proton, then you have just built a hydrogen atom...

As the simulation progresses, you can take data on how
a) Force between the two charges varies with distance between charges;
b) Potential energy changes with distance between charges;
c) Force depends on permittivity.

In each case, take 8 to 10 data points.  Plot your results by hand or by any plotting program.

## EXTENDING THE MODEL

Assign a fixed position to the proton (q1), i.e., make it independent of the mouse position. Assign a variable to its magnitude.

Now create another charge of the breed "centers", and assign a fixed position to it in the graphics window.  Run the model for different positions, magnitude and signs (i.e., "+"ve or "-"ve) of the new "center".

Create many test-charges.  Then place the two "centers", of opposite signs and comparable magnitudes, near the two horizontal edges of the world.  Now run the model.

## RELATED MODELS

* Gravitation

## NETLOGO FEATURES

When a particle moves off of the edge of the world, it doesn't re-appear by wrapping onto the other side (as in most other NetLogo models). The model stops when the particle exits the world.

## CREDITS AND REFERENCES

This model is a part of the NIELS curriculum. The NIELS curriculum has been and is currently under development at Northwestern's Center for Connected Learning and Computer-Based Modeling and the Mind, Matter and Media Lab at Vanderbilt University. For more information about the NIELS curriculum please refer to http://ccl.northwestern.edu/NIELS/.

## HOW TO CITE

If you mention this model or the NetLogo software in a publication, we ask that you include the citations below.

For the model itself:

* Sengupta, P. and Wilensky, U. (2005).  NetLogo Electrostatics model.  http://ccl.northwestern.edu/netlogo/models/Electrostatics.  Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

Please cite the NetLogo software as:

* Wilensky, U. (1999). NetLogo. http://ccl.northwestern.edu/netlogo/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

To cite the NIELS curriculum as a whole, please use:

* Sengupta, P. and Wilensky, U. (2008). NetLogo NIELS curriculum. http://ccl.northwestern.edu/NIELS/. Center for Connected Learning and Computer-Based Modeling, Northwestern University, Evanston, IL.

## COPYRIGHT AND LICENSE

Copyright 2005 Pratim Sengupta and Uri Wilensky.

![CC BY-NC-SA 3.0](http://ccl.northwestern.edu/images/creativecommons/byncsa.png)

This work is licensed under the Creative Commons Attribution-NonCommercial-ShareAlike 3.0 License.  To view a copy of this license, visit https://creativecommons.org/licenses/by-nc-sa/3.0/ or send a letter to Creative Commons, 559 Nathan Abbott Way, Stanford, California 94305, USA.

Commercial licenses are also available. To inquire about commercial licenses, please contact Uri Wilensky at uri@northwestern.edu.

To use this model for academic or commercial research, please contact Pratim Sengupta at <pratim.sengupta@vanderbilt.edu> or Uri Wilensky at <uri@northwestern.edu> for a mutual agreement prior to usage.

<!-- 2005 NIELS Cite: Sengupta, P. -->
@#$#@#$#@
default
true
0
Polygon -7500403 true true 150 5 40 250 150 205 260 250

airplane
true
0
Polygon -7500403 true true 150 0 135 15 120 60 120 105 15 165 15 195 120 180 135 240 105 270 120 285 150 270 180 285 210 270 165 240 180 180 285 195 285 165 180 105 180 60 165 15

arrow
true
0
Polygon -7500403 true true 150 0 0 150 105 150 105 293 195 293 195 150 300 150

box
false
0
Polygon -7500403 true true 150 285 285 225 285 75 150 135
Polygon -7500403 true true 150 135 15 75 150 15 285 75
Polygon -7500403 true true 15 75 15 225 150 285 150 135
Line -16777216 false 150 285 150 135
Line -16777216 false 150 135 15 75
Line -16777216 false 150 135 285 75

bug
true
0
Circle -7500403 true true 96 182 108
Circle -7500403 true true 110 127 80
Circle -7500403 true true 110 75 80
Line -7500403 true 150 100 80 30
Line -7500403 true 150 100 220 30

butterfly
true
0
Polygon -7500403 true true 150 165 209 199 225 225 225 255 195 270 165 255 150 240
Polygon -7500403 true true 150 165 89 198 75 225 75 255 105 270 135 255 150 240
Polygon -7500403 true true 139 148 100 105 55 90 25 90 10 105 10 135 25 180 40 195 85 194 139 163
Polygon -7500403 true true 162 150 200 105 245 90 275 90 290 105 290 135 275 180 260 195 215 195 162 165
Polygon -16777216 true false 150 255 135 225 120 150 135 120 150 105 165 120 180 150 165 225
Circle -16777216 true false 135 90 30
Line -16777216 false 150 105 195 60
Line -16777216 false 150 105 105 60

car
false
0
Polygon -7500403 true true 300 180 279 164 261 144 240 135 226 132 213 106 203 84 185 63 159 50 135 50 75 60 0 150 0 165 0 225 300 225 300 180
Circle -16777216 true false 180 180 90
Circle -16777216 true false 30 180 90
Polygon -16777216 true false 162 80 132 78 134 135 209 135 194 105 189 96 180 89
Circle -7500403 true true 47 195 58
Circle -7500403 true true 195 195 58

circle
false
0
Circle -7500403 true true 0 0 300

circle 2
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240

cow
false
0
Polygon -7500403 true true 200 193 197 249 179 249 177 196 166 187 140 189 93 191 78 179 72 211 49 209 48 181 37 149 25 120 25 89 45 72 103 84 179 75 198 76 252 64 272 81 293 103 285 121 255 121 242 118 224 167
Polygon -7500403 true true 73 210 86 251 62 249 48 208
Polygon -7500403 true true 25 114 16 195 9 204 23 213 25 200 39 123

cylinder
false
0
Circle -7500403 true true 0 0 300

dot
false
0
Circle -7500403 true true 90 90 120

face happy
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 255 90 239 62 213 47 191 67 179 90 203 109 218 150 225 192 218 210 203 227 181 251 194 236 217 212 240

face neutral
false
0
Circle -7500403 true true 8 7 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Rectangle -16777216 true false 60 195 240 225

face sad
false
0
Circle -7500403 true true 8 8 285
Circle -16777216 true false 60 75 60
Circle -16777216 true false 180 75 60
Polygon -16777216 true false 150 168 90 184 62 210 47 232 67 244 90 220 109 205 150 198 192 205 210 220 227 242 251 229 236 206 212 183

fish
false
0
Polygon -1 true false 44 131 21 87 15 86 0 120 15 150 0 180 13 214 20 212 45 166
Polygon -1 true false 135 195 119 235 95 218 76 210 46 204 60 165
Polygon -1 true false 75 45 83 77 71 103 86 114 166 78 135 60
Polygon -7500403 true true 30 136 151 77 226 81 280 119 292 146 292 160 287 170 270 195 195 210 151 212 30 166
Circle -16777216 true false 215 106 30

flag
false
0
Rectangle -7500403 true true 60 15 75 300
Polygon -7500403 true true 90 150 270 90 90 30
Line -7500403 true 75 135 90 135
Line -7500403 true 75 45 90 45

flower
false
0
Polygon -10899396 true false 135 120 165 165 180 210 180 240 150 300 165 300 195 240 195 195 165 135
Circle -7500403 true true 85 132 38
Circle -7500403 true true 130 147 38
Circle -7500403 true true 192 85 38
Circle -7500403 true true 85 40 38
Circle -7500403 true true 177 40 38
Circle -7500403 true true 177 132 38
Circle -7500403 true true 70 85 38
Circle -7500403 true true 130 25 38
Circle -7500403 true true 96 51 108
Circle -16777216 true false 113 68 74
Polygon -10899396 true false 189 233 219 188 249 173 279 188 234 218
Polygon -10899396 true false 180 255 150 210 105 210 75 240 135 240

house
false
0
Rectangle -7500403 true true 45 120 255 285
Rectangle -16777216 true false 120 210 180 285
Polygon -7500403 true true 15 120 150 15 285 120
Line -16777216 false 30 120 270 120

leaf
false
0
Polygon -7500403 true true 150 210 135 195 120 210 60 210 30 195 60 180 60 165 15 135 30 120 15 105 40 104 45 90 60 90 90 105 105 120 120 120 105 60 120 60 135 30 150 15 165 30 180 60 195 60 180 120 195 120 210 105 240 90 255 90 263 104 285 105 270 120 285 135 240 165 240 180 270 195 240 210 180 210 165 195
Polygon -7500403 true true 135 195 135 240 120 255 105 255 105 285 135 285 165 240 165 195

line
true
0
Line -7500403 true 150 0 150 300

line half
true
0
Line -7500403 true 150 0 150 150

pentagon
false
0
Polygon -7500403 true true 150 15 15 120 60 285 240 285 285 120

person
false
0
Circle -7500403 true true 110 5 80
Polygon -7500403 true true 105 90 120 195 90 285 105 300 135 300 150 225 165 300 195 300 210 285 180 195 195 90
Rectangle -7500403 true true 127 79 172 94
Polygon -7500403 true true 195 90 240 150 225 180 165 105
Polygon -7500403 true true 105 90 60 150 75 180 135 105

plant
false
0
Rectangle -7500403 true true 135 90 165 300
Polygon -7500403 true true 135 255 90 210 45 195 75 255 135 285
Polygon -7500403 true true 165 255 210 210 255 195 225 255 165 285
Polygon -7500403 true true 135 180 90 135 45 120 75 180 135 210
Polygon -7500403 true true 165 180 165 210 225 180 255 120 210 135
Polygon -7500403 true true 135 105 90 60 45 45 75 105 135 135
Polygon -7500403 true true 165 105 165 135 225 105 255 45 210 60
Polygon -7500403 true true 135 90 120 45 150 15 180 45 165 90

square
false
0
Rectangle -7500403 true true 30 30 270 270

square 2
false
0
Rectangle -7500403 true true 30 30 270 270
Rectangle -16777216 true false 60 60 240 240

star
false
0
Polygon -7500403 true true 151 1 185 108 298 108 207 175 242 282 151 216 59 282 94 175 3 108 116 108

target
false
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 30 30 240
Circle -7500403 true true 60 60 180
Circle -16777216 true false 90 90 120
Circle -7500403 true true 120 120 60

tree
false
0
Circle -7500403 true true 118 3 94
Rectangle -6459832 true false 120 195 180 300
Circle -7500403 true true 65 21 108
Circle -7500403 true true 116 41 127
Circle -7500403 true true 45 90 120
Circle -7500403 true true 104 74 152

triangle
false
0
Polygon -7500403 true true 150 30 15 255 285 255

triangle 2
false
0
Polygon -7500403 true true 150 30 15 255 285 255
Polygon -16777216 true false 151 99 225 223 75 224

truck
false
0
Rectangle -7500403 true true 4 45 195 187
Polygon -7500403 true true 296 193 296 150 259 134 244 104 208 104 207 194
Rectangle -1 true false 195 60 195 105
Polygon -16777216 true false 238 112 252 141 219 141 218 112
Circle -16777216 true false 234 174 42
Rectangle -7500403 true true 181 185 214 194
Circle -16777216 true false 144 174 42
Circle -16777216 true false 24 174 42
Circle -7500403 false true 24 174 42
Circle -7500403 false true 144 174 42
Circle -7500403 false true 234 174 42

turtle
true
0
Polygon -10899396 true false 215 204 240 233 246 254 228 266 215 252 193 210
Polygon -10899396 true false 195 90 225 75 245 75 260 89 269 108 261 124 240 105 225 105 210 105
Polygon -10899396 true false 105 90 75 75 55 75 40 89 31 108 39 124 60 105 75 105 90 105
Polygon -10899396 true false 132 85 134 64 107 51 108 17 150 2 192 18 192 52 169 65 172 87
Polygon -10899396 true false 85 204 60 233 54 254 72 266 85 252 107 210
Polygon -7500403 true true 119 75 179 75 209 101 224 135 220 225 175 261 128 261 81 224 74 135 88 99

wheel
false
0
Circle -7500403 true true 3 3 294
Circle -16777216 true false 30 30 240
Line -7500403 true 150 285 150 15
Line -7500403 true 15 150 285 150
Circle -7500403 true true 120 120 60
Line -7500403 true 216 40 79 269
Line -7500403 true 40 84 269 221
Line -7500403 true 40 216 269 79
Line -7500403 true 84 40 221 269

x
false
0
Polygon -7500403 true true 270 75 225 30 30 225 75 270
Polygon -7500403 true true 30 75 75 30 270 225 225 270
@#$#@#$#@
NetLogo 6.1.1
@#$#@#$#@
need-to-manually-make-preview-for-this-model
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
@#$#@#$#@
default
0.0
-0.2 0 0.0 1.0
0.0 1 1.0 0.0
0.2 0 0.0 1.0
link direction
true
0
Line -7500403 true 150 150 90 180
Line -7500403 true 150 150 210 180
@#$#@#$#@
0
@#$#@#$#@
