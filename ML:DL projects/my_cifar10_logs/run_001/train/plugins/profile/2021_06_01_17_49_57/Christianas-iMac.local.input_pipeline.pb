	y!?x@y!?x@!y!?x@	??%?"?????%?"???!??%?"???"{
=type.googleapis.com/tensorflow.profiler.PerGenericStepDetails:y!?x@2?F? ??A=b???@Y??YKi??rEagerKernelExecute 0*	?(\?¥V@2v
?Iterator::Model::ParallelMapV2::Zip[0]::FlatMap[0]::Concatenate? [??ˠ?!~=??B@)??u?ӝ?1b??G?@@:Preprocessing2f
/Iterator::Model::ParallelMapV2::Zip[0]::FlatMap<?D???!L?ů??N@)?S?*???1?݊ZT}9@:Preprocessing2U
Iterator::Model::ParallelMapV2?7? ???!?o????*@)?7? ???1?o????*@:Preprocessing2l
5Iterator::Model::ParallelMapV2::Zip[1]::ForeverRepeatk?#?]J??!?6B?/@)?d??1??}J??)@:Preprocessing2F
Iterator::ModelOI?V??!??ذ2@)2??|?s?16(u?!@:Preprocessing2?
OIterator::Model::ParallelMapV2::Zip[0]::FlatMap[0]::Concatenate[0]::TensorSlice%u?n?!??2?:@)%u?n?1??2?:@:Preprocessing2Z
#Iterator::Model::ParallelMapV2::Zip??-@۲?!?????ST@)?r?9>Zl?1????g?@:Preprocessing2x
AIterator::Model::ParallelMapV2::Zip[1]::ForeverRepeat::FromTensorNA~6r?d?!?0c?~@)NA~6r?d?1?0c?~@:Preprocessing:?
]Enqueuing data: you may want to combine small input data chunks into fewer but larger chunks.
?Data preprocessing: you may increase num_parallel_calls in <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#map" target="_blank">Dataset map()</a> or preprocess the data OFFLINE.
?Reading data from files in advance: you may tune parameters in the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch size</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave cycle_length</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer_size</a>)
?Reading data from files on demand: you should read data IN ADVANCE using the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer</a>)
?Other data reading or processing: you may consider using the <a href="https://www.tensorflow.org/programmers_guide/datasets" target="_blank">tf.data API</a> (if you are not using it now)?
:type.googleapis.com/tensorflow.profiler.BottleneckAnalysis?
device?Your program is NOT input-bound because only 1.5% of the total step time sampled is waiting for input. Therefore, you should focus on reducing other time.no*no9??%?"???I?hu/?X@Zno#You may skip the rest of this page.B?
@type.googleapis.com/tensorflow.profiler.GenericStepTimeBreakdown?
	2?F? ??2?F? ??!2?F? ??      ??!       "      ??!       *      ??!       2	=b???@=b???@!=b???@:      ??!       B      ??!       J	??YKi????YKi??!??YKi??R      ??!       Z	??YKi????YKi??!??YKi??b      ??!       JCPU_ONLYY??%?"???b q?hu/?X@