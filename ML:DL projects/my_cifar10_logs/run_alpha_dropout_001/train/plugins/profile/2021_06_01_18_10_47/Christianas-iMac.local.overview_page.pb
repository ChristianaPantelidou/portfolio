?	?t?y?>@?t?y?>@!?t?y?>@	mgB?S??mgB?S??!mgB?S??"{
=type.googleapis.com/tensorflow.profiler.PerGenericStepDetails:?t?y?>@=???????A`?+?۶@Y?5C?(??rEagerKernelExecute 0*	??"??VZ@2l
5Iterator::Model::ParallelMapV2::Zip[1]::ForeverRepeat3Mg'???!?qHg;?I@)????c??1)?d?uH@:Preprocessing2U
Iterator::Model::ParallelMapV2??gx???!????<@)??gx???1????<@:Preprocessing2v
?Iterator::Model::ParallelMapV2::Zip[0]::FlatMap[0]::Concatenatez??y???!?ˌ??u#@)??.?.y?1???wUW@:Preprocessing2F
Iterator::Model;?p?GR??!i?I`?@@),??ypwv?1??<D??@:Preprocessing2?
OIterator::Model::ParallelMapV2::Zip[0]::FlatMap[0]::Concatenate[0]::TensorSlice%=?N?p?!6???'@)%=?N?p?16???'@:Preprocessing2Z
#Iterator::Model::ParallelMapV2::Zip䃞ͪϱ?!̯~?O?P@)?v??-up?1/?A<O?@:Preprocessing2x
AIterator::Model::ParallelMapV2::Zip[1]::ForeverRepeat::FromTensorA?
??a?!?,?? @)A?
??a?1?,?? @:Preprocessing2f
/Iterator::Model::ParallelMapV2::Zip[0]::FlatMap ???!6??!IP?o?p&@)2?g@?Y?1#??????:Preprocessing:?
]Enqueuing data: you may want to combine small input data chunks into fewer but larger chunks.
?Data preprocessing: you may increase num_parallel_calls in <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#map" target="_blank">Dataset map()</a> or preprocess the data OFFLINE.
?Reading data from files in advance: you may tune parameters in the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch size</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave cycle_length</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer_size</a>)
?Reading data from files on demand: you should read data IN ADVANCE using the following tf.data API (<a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#prefetch" target="_blank">prefetch</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/Dataset#interleave" target="_blank">interleave</a>, <a href="https://www.tensorflow.org/api_docs/python/tf/data/TFRecordDataset#class_tfrecorddataset" target="_blank">reader buffer</a>)
?Other data reading or processing: you may consider using the <a href="https://www.tensorflow.org/programmers_guide/datasets" target="_blank">tf.data API</a> (if you are not using it now)?
:type.googleapis.com/tensorflow.profiler.BottleneckAnalysis?
device?Your program is NOT input-bound because only 1.6% of the total step time sampled is waiting for input. Therefore, you should focus on reducing other time.no*no9mgB?S??Ib?????X@Zno#You may skip the rest of this page.B?
@type.googleapis.com/tensorflow.profiler.GenericStepTimeBreakdown?
	=???????=???????!=???????      ??!       "      ??!       *      ??!       2	`?+?۶@`?+?۶@!`?+?۶@:      ??!       B      ??!       J	?5C?(???5C?(??!?5C?(??R      ??!       Z	?5C?(???5C?(??!?5C?(??b      ??!       JCPU_ONLYYmgB?S??b qb?????X@Y      Y@qT?芙C1@"?	
device?Your program is NOT input-bound because only 1.6% of the total step time sampled is waiting for input. Therefore, you should focus on reducing other time.b
`input_pipeline_analyzer (especially Section 3 for the breakdown of input operations on the Host)Q
Otf_data_bottleneck_analysis (find the bottleneck in the tf.data input pipeline)m
ktrace_viewer (look at the activities on the timeline of each Host Thread near the bottom of the trace view)"T
Rtensorflow_stats (identify the time-consuming operations executed on the CPU_ONLY)"Z
Xtrace_viewer (look at the activities on the timeline of each CPU_ONLY in the trace view)*?
?<a href="https://www.tensorflow.org/guide/data_performance_analysis" target="_blank">Analyze tf.data performance with the TF Profiler</a>*y
w<a href="https://www.tensorflow.org/guide/data_performance" target="_blank">Better performance with the tf.data API</a>2M
=type.googleapis.com/tensorflow.profiler.GenericRecommendation
nono2no:
Refer to the TF2 Profiler FAQb?17.3% of Op time on the host used eager execution. Performance could be improved with <a href="https://www.tensorflow.org/guide/function" target="_blank">tf.function.</a>2"CPU: B 