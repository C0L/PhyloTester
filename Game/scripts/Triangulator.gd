extends Spatial
var startPos
var main

var TriObj0 = preload("../scenes/TriObj0.tscn")
var TriObj1 = preload("../scenes/TriObj1.tscn")
var TriObj2 = preload("../scenes/TriObj2.tscn")

var bob = Vector3(0, .001, 0)

# Called when the node enters the scene tree for the first time.
func _ready():
	startPos = self.get_translation()
	$StaticWeapons/w0/CSGMesh0.show()
	main = $StaticWeapons/w0/CSGMesh0



func setMain(num):
	#Hide all
	$StaticWeapons/w0/CSGMesh0.hide()
	$StaticWeapons/w0/CSGMesh1.hide()
	$StaticWeapons/w0/CSGMesh2.hide()
	
	if num == 0:
		main = $StaticWeapons/w0/CSGMesh0
	elif num == 1:
		main = $StaticWeapons/w0/CSGMesh1
	elif num == 2:
		main = $StaticWeapons/w0/CSGMesh2
	main.show()

func process_input():
	if Input.is_action_pressed("weap0"):
		setMain(0)
	elif Input.is_action_pressed("weap1"):
		setMain(1)
	elif Input.is_action_pressed("weap2"):
		setMain(2)

	if Input.is_action_just_pressed("lmbutton"):
		main.hide()
		var clone
		if main == $StaticWeapons/w0/CSGMesh0:
			clone = TriObj0.instance()
		elif main == $StaticWeapons/w0/CSGMesh1:
			clone = TriObj1.instance()
		else:
			clone = TriObj2.instance()
	
		var scene_root = get_tree().root.get_children()[0]
		scene_root.add_child(clone)
		clone.global_transform = self.global_transform

#func triangulate():
#	var clone = $TriOb.instance()
#	var scene_root = get_tree().root.get_children()[0]
#	scene_root.add_child(clone)
#	clone.global

# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta):
	process_input()
	
	var nPos = self.get_translation()
	if (nPos[1] > (startPos[1] + .3)):
		bob = Vector3(0, -.001, 0)
	elif (nPos[1] < (startPos[1] - .3)):
		bob = Vector3(0, .001, 0)
	
	main.rotate_x(-.005)
	main.rotate_y(.01)
	main.rotate_z(-.005)
	self.translate(bob)
