extends Spatial

var FIRE_SPEED = 35
#var DAMAGE = 15

const KILL_TIMER = 50
var timer = 0

var hit_something = false

func _ready():
	#pass
	$Area.connect("body_entered", self, "collided")


func _physics_process(delta):
	var forward_dir = -global_transform.basis.z.normalized()
	global_translate(forward_dir * FIRE_SPEED * delta)

	timer += delta
	if timer >= KILL_TIMER:
		queue_free()


func collided(body):
	# Ignore if touching the static weapons
	if body.get_name() == "StaticWeapons":
		return

	#if hit_something == false:
		#if body.has_method("bullet_hit"):
			#ody.bullet_hit(BULLET_DAMAGE, global_transform)

	#hit_something = true

	queue_free()
	

# Called every frame. 'delta' is the elapsed time since the previous frame.
func _process(delta):
	#var nPos = self.get_translation()
	#if (nPos[1] > (startPos[1] + .3)):
	#	bob = Vector3(0, -.001, 0)
	#elif (nPos[1] < (startPos[1] - .3)):
	#	bob = Vector3(0, .001, 0)
	
	$CSGMesh.rotate_x(-.005)
	$CSGMesh.rotate_y(.01)
	$CSGMesh.rotate_z(-.005)
	#self.translate(bob)
