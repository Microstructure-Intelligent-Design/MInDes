#pragma once
#ifdef _WIN32
#include <Windows.h>
#include <process.h>
#include <psapi.h>
#else
#include <unistd.h>
#include <sys/types.h>
#include <unistd.h>  // getpid(), fork(), etc.
#endif
#include "../input_modules/ioFiles_Params.h"
#include "../postprocess_modules/ShowLoopInfo.h"
#include "../Module.h"
#include <cstring>
namespace pf {
	namespace cpu_memory_usage {
        inline int64_t last_time = 0;
        inline int64_t last_system_time = 0;
        // get current process pid
        int GetCurrentPid();
        // get specific process cpu occupation ratio by pid
// FIXME: can also get cpu and mem status from popen cmd
// the info line num in /proc/{pid}/status file
#define VMRSS_LINE 22
#define PROCESS_ITEM 14
#ifdef WIN32
        uint64_t convert_time_format(const FILETIME* ftime);
#else
        const char* get_items(const char* buffer, unsigned int item);
        unsigned long get_cpu_total_occupy();
        unsigned long get_cpu_proc_occupy(int pid);
#endif
        float GetCpuUsageRatio(int pid);
        double GetMemoryUsage(int pid);
        void exec_pre_iii();
        void exec_post_ii();
        void deinit();

        void init_cpu_memory_usage();
	}
}