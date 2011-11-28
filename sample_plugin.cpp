// created by Sebastian Reiter
// s.b.reiter@googlemail.com
// 12.09.2011 (m,d,y)

#include "registry/registry.h"
#include "common/ug_config.h"
#include "common/error.h"
#include <string>

using namespace std;

#ifdef UG_PARALLEL
#include "pcl/pcl_util.h"
#endif

void PluginSaysHello()
{
#ifdef UG_PARALLEL
	pcl::SynchronizeProcesses();
	cout << "Hello, I'm your plugin on proc " <<
				pcl::GetProcRank() << "." << endl;
	pcl::SynchronizeProcesses();
#else
	UG_LOG("Hello, I'm your personal plugin in serial environment!\n");
#endif
}

void CrashFct(string reason)
{
	UG_THROW("I Crash because: "<< reason);
}

void CrashFctFatal(string reason)
{
	UG_THROW_FATAL("I Crash fatal because: "<< reason);
}

void PluginCrashes()
{
	try{
		CrashFct("Some funny reason");
	}
	catch(bad_alloc err)
	{
		UG_LOG("Bad Alloc");
	}
	UG_CATCH_THROW("Something wrong in PluginCrashes");
}

void PluginCrashesFatal()
{
	try{
		CrashFctFatal("Some fatal reason");
	}
	catch(bad_alloc err)
	{
		UG_LOG("Bad Alloc");
	}
	UG_CATCH_THROW("Something wrong in PluginCrashesFatal");
}

extern "C" UG_API void InitUGPlugin(ug::bridge::Registry* reg, string parentGroup)
{
	string grp(parentGroup); grp.append("SamplePlugin/");

	reg->add_function("PluginSaysHello", &PluginSaysHello, grp)
		.add_function("PluginCrashes", &PluginCrashes, grp)
		.add_function("PluginCrashesFatal", &PluginCrashesFatal, grp);
}
