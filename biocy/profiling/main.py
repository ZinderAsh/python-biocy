import cProfile
import profile_speed
import profile_memory

cProfile.run('profile_speed.main()')
profile_memory.main()

