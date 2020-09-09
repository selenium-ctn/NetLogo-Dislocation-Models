breed [atoms atom]
breed [fl-ends fl-end] ; turtles at the ends of the force lines, show which direction the force is acting
undirected-link-breed [fl-links fl-link] ; force line links
undirected-link-breed [atom-links atom-link] ; links between atoms

atoms-own [
  fx     ; x-component of force vector
  fy     ; y-component of force vector
  vx     ; x-component of velocity vector
  vy     ; y-component of velocity vector
  mass   ; mass of atom
  pinned? ; False if the atom isn't pinned in place, True if it is (for boundaries)
  ex-force-applied? ; is an external force directly applied to this atom? False if no, True if yes
  total-PE
  my-new-fx
  my-new-fy
]

globals [
  eps ; used in LJ force. Well depth; measure of how strongly particles attract each other
  sigma ; used in LJ force. Distance at which intermolecular potential between 2 particles is 0
  cutoff-dist ; each atom is influenced by its neighbors within this distance (LJ force)
  dt ; time step for the velocity verlet algorithm
  sqrt-2-kb-over-m  ; constant. Used when calculating the thermal velocity. Square root of (2 * boltzmann constant / m). It is arbitrary for this simulation since the units are also arbitrary.
  link-check-dist ; each atom links with neighbors within this distance
  prev-lattice-view ; the lattice view in the previous time step
  upper-left-fl ; upper left force line - shear
  left-fl ; left force line - tension, compression
  right-edge ; where the right side of the sample is (xcor) - tension, compression - used in determining length of sample
  orig-length ; original length of sample
  prev-length ; length of sample in previous time step
  median-ycor ; median ycor of atoms from initial lattice setup
  top-neck-atoms ; agentset of atoms on the top of the neck (thin region) (tension). Used in calculating stress
  bottom-neck-atoms ; agentset of atoms on the bottom of the neck (thin region) (tension). Used in calculating stress
  num-forced-atoms ; number of atoms receiving external force directly
  unpinned-atoms ; atoms that are not pinned
  total-ex-force
  indiv-force
]

;;;;;;;;;;;;;;;;;;;;;;
;; Setup Procedures ;;
;;;;;;;;;;;;;;;;;;;;;;

to setup
  clear-all
  set eps .07
  set sigma .907
  set cutoff-dist 5
  set dt .1
  set sqrt-2-kb-over-m (1 / 50)
  set link-check-dist 1.5
  setup-atoms-and-links-and-fls
  init-velocity
  update-lattice-view
  reset-ticks
end

to setup-atoms-and-links-and-fls
  if force-mode = "Tension" and atoms-per-column mod 2 = 0 [ set atoms-per-column atoms-per-column + 1] ; making a symmetrical sample for tension mode
  create-atoms atoms-per-row * atoms-per-column [
    set shape "circle"
    set color blue
    set mass 1
    set pinned? False
    set ex-force-applied? False
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
      set pinned? True
    ]
    ]
    force-mode = "Tension"[
      ask atoms with [xcor = xmin] [die] ; creating the symmetrical shape
      set xmin min [xcor] of atoms
      ask atoms with [(ycor >= ymax - 1 or ycor <= ymin + 1) and xcor <= xmax - 3.5 and xcor >= xmin + 3.5] [die]
      ask atoms with [xcor = xmax or xcor = xmax - .5 ] [set pinned? True]

      set top-neck-atoms atoms with [xcor <= xmax - 3.5 and xcor >= xmin + 3.5] with-max [ycor] ; defining top and bottom neck agentsets
      set bottom-neck-atoms atoms with [xcor <= xmax - 3.5 and xcor >= xmin + 3.5] with-min [ycor]
      ask atoms with [ xcor >= xmin and xcor <= xmin + 3 ][
        set ex-force-applied? True
        set shape "circle-dot"
      ]
      set num-forced-atoms count atoms with [ex-force-applied?]
    ]
    force-mode = "Compression" [
      ask atoms with [xcor = xmax or xcor = xmax - .5 ] [set pinned? True]
    ]
  )

  if create-dislocation? [ ; creating the dislocation
    let curr-y-cor median-ycor
    let curr-x-cor median-xcor
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
  ask atom-links [ ; stylizing/coloring links
    color-links
  ]

  (ifelse force-mode = "Tension"  [ ; set up force lines
    create-fl-ends 2
    set left-fl xmin
    set right-edge xmax
    set orig-length right-edge - left-fl
    ask one-of fl-ends with [xcor = 0 and ycor = 0] [
      set xcor left-fl
      set ycor ymax + 2 ]
    ask one-of fl-ends with [xcor = 0 and ycor = 0] [
      set xcor left-fl
      set ycor ymin - 2
      create-fl-link-with one-of other fl-ends with [xcor = left-fl]]
    ask fl-ends [
      set color white
      set heading 270
    ]
    if force-mode = "Tension" [
      set prev-length orig-length
    ]
  ]
    force-mode = "Compression" [
      create-fl-ends 2
      set left-fl xmin
      set right-edge xmax
      set orig-length right-edge - left-fl
      ask one-of fl-ends with [xcor = 0 and ycor = 0] [
        set xcor left-fl
        set ycor max-pycor - 2 ]
      ask one-of fl-ends with [xcor = 0 and ycor = 0] [
        set xcor left-fl
        set ycor min-pycor + 2
        create-fl-link-with one-of other fl-ends with [xcor = left-fl]]
      ask fl-ends [
        set color white
        set heading 90
      ]
    ]
    force-mode = "Shear" [
      create-fl-ends 2
      set upper-left-fl min [xcor] of atoms with [ ycor >= median-ycor ]
      ask one-of fl-ends with [xcor = 0 and ycor = 0] [
        set xcor upper-left-fl
        set ycor ymax + 2 ]
      ask one-of fl-ends with [xcor = 0 and ycor = 0] [
        set xcor upper-left-fl
        set ycor median-ycor
        hide-turtle
        create-fl-link-with one-of other fl-ends]
      ask fl-ends [
        set color white
        set heading 90 ]
  ])
  ask fl-links [
    set color white
  ]
  ask atoms with [pinned?] [ set shape "circle-x"]
  set unpinned-atoms atoms with [not pinned?]
end

to init-velocity ; initializes velocity for each atom based on the initial system-temp. Creates a random aspect in the
                 ; velocity split between the x velocity and the y velocity
  let speed-avg sqrt-2-kb-over-m * sqrt system-temp
  ask unpinned-atoms [
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
  set total-ex-force 0
  if lattice-view != prev-lattice-view [ update-lattice-view ]
  control-temp
  ask atom-links [die]
  if mouse-down? [ delete-atoms ]
  ask unpinned-atoms [ ; moving happens before velocity and force update in accordance with velocity verlet
    move
  ]
  calculate-fl-positions
  if force-mode = "Tension" and auto-increment-force? [ adjust-force ]
  identify-force-atoms
  ask atoms [
    update-force-and-velocity-and-links
  ]
  ifelse f-app - total-ex-force > 0 [
    set indiv-force ( f-app - total-ex-force ) / num-forced-atoms ]
  [set indiv-force 0 ]
  ask atoms [update-2]
  ask atom-links [ ; stylizing/coloring links
    color-links
  ]
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
  set prev-lattice-view lattice-view
end

to control-temp ; this heats or cools the system based on the average temperature of the system compared to the set system-temp
  let current-speed-avg mean [ sqrt (vx ^ 2 + vy ^ 2) ] of unpinned-atoms
  let target-speed-avg sqrt-2-kb-over-m * sqrt system-temp
  let scaling-factor target-speed-avg / current-speed-avg
  if current-speed-avg != 0 [
    ask unpinned-atoms [
      set vx vx * scaling-factor
      set vy vy * scaling-factor
    ]
  ]
end



to move  ; atom procedure, uses velocity-verlet algorithm
  set xcor velocity-verlet-pos xcor vx (fx / mass)
  set ycor velocity-verlet-pos ycor vy (fy / mass)
  if xcor > max-pxcor or xcor < min-pxcor [
    die ; kills atoms when they move off the world
  ]
end

to calculate-fl-positions ; (calculate new force line positions)
  ifelse force-mode = "Shear" [
    set upper-left-fl min [xcor] of atoms with [ ycor >= median-ycor ]
    ask fl-ends [ set xcor upper-left-fl]
  ]
  [ ; force-mode = tension or compression
    set left-fl min [xcor] of atoms
    ask fl-ends with [xcor < 0] [ set xcor left-fl]
  ]
end

to identify-force-atoms ; (find the atoms closest to the force line that will be the ones receiving the external force)
  (ifelse force-mode = "Shear" [
    ask atoms [ set ex-force-applied?  False ]
    let forced-atoms atoms with [ ycor >= median-ycor and (distancexy upper-left-fl ycor) <= 1]
    set num-forced-atoms count forced-atoms
    ask forced-atoms [
      set ex-force-applied?  True
    ]
    ]
    force-mode = "Compression" [
      ask atoms [ set ex-force-applied?  False ]
      let forced-atoms atoms with [ (distancexy left-fl ycor) <= 1]
      set num-forced-atoms count forced-atoms
      ask forced-atoms [
        set ex-force-applied?  True
    ]
  ]) ; for tension, the same atoms in the left shoulder of the sample always receive the force
end


to adjust-force
  if precision prev-length 6 >= precision (right-edge - left-fl) 6 [ set f-app precision (f-app + .005) 3 ]
  ; increments f-app-auto if the sample has reached an equilibrium or if the previous sample length is greater than the current sample length
  set prev-length (right-edge - left-fl)
end

to update-force-and-velocity-and-links
  let new-fx 0
  let new-fy 0
  let total-potential-energy 0
  let in-radius-atoms other atoms in-radius cutoff-dist
  ask in-radius-atoms [ ; each atom calculates the force it feels from its neighboring atoms and sums these forces
    let r distance myself
    let indiv-PE-and-force (LJ-poten-and-force r)
    let force item 1 indiv-PE-and-force
    set total-potential-energy total-potential-energy + item 0 indiv-PE-and-force
    face myself
    rt 180
    set new-fx new-fx + (force * dx)
    set new-fy new-fy + (force * dy)
    ]
  set total-PE total-potential-energy

  set my-new-fx new-fx
  set my-new-fy new-fy

  if ex-force-applied? and my-new-fx > 0 [ set total-ex-force total-ex-force + new-fx
  set my-new-fx 0 ]

end

to update-2
  if not pinned? [
    ; adjusting the forces to account for any external applied forces
    let ex-force 0
    if ex-force-applied? [ set ex-force report-new-force
    set my-new-fy 0]
    if shape = "circle-dot" and not ex-force-applied? [ set shape "circle" ]
    set my-new-fx ex-force + my-new-fx

    ; updating velocity and force
    set vx velocity-verlet-velocity vx (fx / mass) (my-new-fx / mass)
    set vy velocity-verlet-velocity vy (fy / mass) (my-new-fy / mass)
    set fx my-new-fx
    set fy my-new-fy
  ]

  update-atom-color total-PE
  ;update-links in-radius-atoms
end

to update-atom-color [total-force] ; updating atom color
  (ifelse update-atom-color? [
    set-color total-force
  ]
   [ set color blue ])
end

to update-links [in-radius-atoms] ; updating links
  if show-diagonal-right-links? [
    set heading 330
    link-with-atoms-in-cone in-radius-atoms
  ]
  if show-diagonal-left-links? [
    set heading 30
    link-with-atoms-in-cone in-radius-atoms
  ]
  if show-horizontal-links? [
    set heading 90
    link-with-atoms-in-cone in-radius-atoms
  ]
end

to link-with-atoms-in-cone [atom-set]
  let in-cone-atoms (atom-set in-cone link-check-dist 60)
    if any? in-cone-atoms [
      create-atom-link-with min-one-of in-cone-atoms [distance myself]
    ]
end

to-report report-new-force
  set shape "circle-dot"
  (ifelse force-mode = "Tension" [
    ;report -1 * f-app / num-forced-atoms
    ;report -1 * indiv-force
    ]
    [ ; Shear and Compression
      report f-app / num-forced-atoms
    ]
  )
end

to-report LJ-poten-and-force [ r ] ; for the force, positive = attractive, negative = repulsive
  let third-pwr (sigma / r) ^ 3
  let sixth-pwr third-pwr ^ 2
  let twelfth-pwr sixth-pwr ^ 2
  let force (-48 * eps / r ) * (twelfth-pwr - (1 / 2) * sixth-pwr) + .0001
  ;report (-48 * eps / r )* ((sigma / r) ^ 12 - (1 / 2) * (sigma / r) ^ 6) + .0001
  let potential (4 * eps * (twelfth-pwr - sixth-pwr)) + .00001
  report list potential force
end

to-report velocity-verlet-pos [pos v a]  ; position, velocity and acceleration
  report pos + v * dt + (1 / 2) * a * (dt ^ 2)
end

to-report velocity-verlet-velocity [v a new-a]  ; velocity, acceleration, new acceleration
  report v + (1 / 2) * (new-a + a) * dt
end

to set-color [v]
  set color scale-color blue v -.9 0
end

to-report strain ; tension only
  report ((right-edge - left-fl) - orig-length) / orig-length
end

to-report stress ; tension only
  let avg-max mean [ycor] of top-neck-atoms
  let avg-min mean [ycor] of bottom-neck-atoms
  let min-A avg-max - avg-min
;  ifelse ( f-app - total-ex-force ) > 0
;  [report 0]
;  [report total-ex-force / min-A ]
  ifelse indiv-force != 0 [
    report (f-app + total-ex-force) / min-A ]
  [report total-ex-force / min-A ]
end

to-report report-indiv-ex-force
  report f-app / num-forced-atoms
end

to-report rep-force
    ifelse indiv-force != 0 [
    report (f-app + total-ex-force) ]
  [report total-ex-force ]
end

to color-links
  set thickness .25 ; necessary bc the links die and reform every tick
  let min-eq-bond-len .995
  let max-eq-bond-len 1.018073
  (ifelse
    link-length < min-eq-bond-len [
      let tmp-len sqrt(min-eq-bond-len - link-length)
      let tmp-color extract-rgb scale-color red tmp-len 1 -.2
      set color insert-item 3 tmp-color (125 + (1 + tmp-len) * 30) ]
    link-length > max-eq-bond-len [
      let tmp-len sqrt (link-length - max-eq-bond-len)
      let tmp-color extract-rgb scale-color yellow tmp-len 1 -.2
      set color insert-item 3 tmp-color (125 + (1 + tmp-len) * 30)]
    [ let tmp-color extract-rgb white
      set color insert-item 3 tmp-color 125 ])
end
@#$#@#$#@
GRAPHICS-WINDOW
197
10
865
679
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
-16
16
-16
16
1
1
1
ticks
30.0

BUTTON
9
257
95
290
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
103
257
188
290
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
1

SLIDER
10
351
182
384
system-temp
system-temp
0
.4
0.296
.001
1
NIL
HORIZONTAL

SLIDER
10
394
182
427
f-app
f-app
0
30
<<<<<<< Updated upstream
1.91
=======
1.005
>>>>>>> Stashed changes
.1
1
N
HORIZONTAL

SWITCH
17
128
176
161
create-dislocation?
create-dislocation?
1
1
-1000

SWITCH
879
13
1043
46
update-atom-color?
update-atom-color?
1
1
-1000

CHOOSER
12
299
179
344
lattice-view
lattice-view
"large-atoms" "small-atoms" "hide-atoms"
1

SWITCH
882
55
1080
88
show-diagonal-right-links?
show-diagonal-right-links?
0
1
-1000

SWITCH
882
95
1074
128
show-diagonal-left-links?
show-diagonal-left-links?
0
1
-1000

SWITCH
881
133
1058
166
show-horizontal-links?
show-horizontal-links?
0
1
-1000

SLIDER
11
169
183
202
atoms-per-row
atoms-per-row
5
20
18.0
1
1
NIL
HORIZONTAL

SLIDER
11
211
183
244
atoms-per-column
atoms-per-column
5
20
13.0
1
1
NIL
HORIZONTAL

MONITOR
881
184
1085
229
external force per forced atom (N)
report-indiv-ex-force
3
1
11

PLOT
881
247
1081
397
Stress-Strain Curve - Tension
strain
stress
0.0
0.05
0.0
0.1
true
false
"" ""
PENS
"default" 1.0 2 -16777216 true "" "if force-mode = \"Tension\" [ plotxy strain stress ]"

MONITOR
884
416
992
461
current f-app (N)
;f-app\nrep-force
5
1
11

MONITOR
997
415
1113
460
sample length (rm)
right-edge - left-fl
5
1
11

TEXTBOX
1125
420
1245
528
<- in terms of the equilibrium interatomic distance between two atoms (rm)
11
0.0
1

TEXTBOX
881
468
1031
486
NIL
11
0.0
1

TEXTBOX
887
467
1037
495
Color Key\nLinks: 
11
0.0
1

TEXTBOX
887
497
1037
515
high compression: dark red
11
13.0
1

TEXTBOX
887
512
1094
530
low compression: light red (+ grey tone)
11
18.0
1

TEXTBOX
886
526
1036
544
equilibrium: grey
11
5.0
1

TEXTBOX
886
539
1089
567
low tension: light yellow (+ grey tone)
11
0.0
1

TEXTBOX
1074
535
1124
555
■■■■
16
48.0
1

TEXTBOX
887
555
1037
573
high tension: dark yellow
11
44.0
1

TEXTBOX
889
579
1039
597
Atoms:
11
0.0
1

TEXTBOX
889
592
1039
610
low potential energy: dark blue 
11
103.0
1

TEXTBOX
889
606
1100
634
high potential energy: light blue (-> white)
11
107.0
1

TEXTBOX
889
618
1085
702
pinned atoms (do not move): black cross\natoms affected by an external force: black dot, near a white line with arrows on the end 
11
0.0
1

SWITCH
13
69
192
102
auto-increment-force?
auto-increment-force?
0
1
-1000

TEXTBOX
50
104
200
122
(^ Tension only)
11
0.0
1

@#$#@#$#@
## WHAT IS IT?

This model allows the user to observe the effects of external forces on a close-packed 2D crystal lattice. It also gives a qualitative image of stress/strain fields around edge dislocations. An edge dislocation can be initialized within the material, and shear, tension, or compression forces can be applied to the material. This is a Molecular Dynamics simulation, meaning an interatomic potential governs the energy and motion of each atom. Here, we are using the Lennard-Jones potential. Atoms attempt to minimize their energy by moving to an equilibrium distance from the atoms surrounding them, which means that the atoms are moved in response to the force they feel from the surrounding atoms. The position and velocity of the atoms are updated every time step with the velocity Verlet algorithm. 

## HOW IT WORKS

The Lennard-Jones potential shows that there is an equilibrium distances between two atoms where the potential energy of each atom is minimized. If the atoms are farther than this equilibrium distance, they will attract each other, and if they are closer, they will repel each other. 

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

circle-dot
true
0
Circle -7500403 true true 0 0 300
Circle -16777216 true false 88 88 124

circle-x
false
0
Circle -7500403 true true 0 0 300
Rectangle -16777216 true false 0 120 315 165
Rectangle -16777216 true false 135 -15 180 300

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
