project(class2)
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_FLAGS "-fPIC")

include_directories("../..")


################################################################ test camera model_my
#set(CAMERA_MODEL_MY
#        class2_test_camera_my.cc
#        )
#
#add_executable(class2_test_camera_my ${CAMERA_MODEL_MY}  )
#target_link_libraries(class2_test_camera_my util core features )


################################################################# test fundamental 8 point_my
#set(FUNDAMENTAL_8_POINT_FILE
#        class2_test_fundamental_8_point_my.cc)
#add_executable(task2_test_fundamental_8_point_my ${FUNDAMENTAL_8_POINT_FILE} )
#target_link_libraries(task2_test_fundamental_8_point_my sfm util core features )




##################################################################### test fundamental ransac_my
set(FUNDAMENTAL_RANSAC_FILE
	       class2_test_fundamental_ransac.my.cc)
add_executable(class2_test_fundamental_ransac.my ${FUNDAMENTAL_RANSAC_FILE})
target_link_libraries(class2_test_fundamental_ransac.my sfm util core features )


##################################################################### test pose from fundamental_my
#set(POSE_FROM_FUNDAMENTAL
#        class2_test_pose_from_fundamental.cc
#        )
#add_executable(task2_test_pose_from_fundamental ${POSE_FROM_FUNDAMENTAL} )
#target_link_libraries(task2_test_pose_from_fundamental sfm util core features )
